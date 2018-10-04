function [ChanMod_output UE_input] = LTE_channel_model(LTE_params,ChanMod,ChanMod_output, BS_output, SNR_vec, UE_input, subframe_i)
% LTE channel model - to filter the output of the transmitter.
% [chan_output] = LTE_channel_model(BS_output, SNR)
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   BS_output           ... [1x1]struct - Base Station output
%                                   [LTE_params.TxSymbols x 1]double transmit signal y_tx,
%                                   [1x1]struct - scheduler output and coding parameters, e.g:
%                                   [1 x # data bits]logical genie.data_bits
%                                   [1 x # sent bits(coded)]logical genie.sent_bits
%           SNR                 ... [1]double current SNR value
% output:   chan_output         ... [1x1]struct - channel output:
%                                   [LTE_params.TxSymbols x 1]double receive signal y_rx
%
% date of creation: 2008/08/11
% last changes: 2008/09/18  Bosanska    changed the name of chan_output -> ChanMod_output
%                                       added input [1 x 1]struct ChanMod with channel model parameters
%                                       added input [1 x 1]struct ChanMod_output with channel matrix H
%               2008/11/04  Bosanska    implemented Oversampling, RRC Transmit Filter, Resampling and Interpolation,
%                                       Resampling and Decimation, RRC Receive Filter, Downsampling
%                                       implemented perfect channel knowledge (Freq Response)
%               2008/11/07  Bosanska    removed Oversampling, RRC Transmit Filter, Resampling and Interpolation,
%                                       Resampling and Decimation, RRC, Receive Filter, Downsampling

%% Channel Model
for uu = 1:LTE_params.nUE*LTE_params.nBS
    
    SNR = SNR_vec(uu);
    for bb = 1:LTE_params.nBS
        switch ChanMod.type
            case {'AWGN','flat Rayleigh'}
                y_tx_filtered = ChanMod_output{bb,uu}.H*BS_output(bb).y_tx.';
                y_tx_filtered = y_tx_filtered.';
            case {'PedA','PedB','PedBcorr','VehA','VehB','TU','RA','HT','EPedA','EVehA','ETU','winner_II','Rayleigh2'}
                switch ChanMod.filtering
                    case 'BlockFading'
                        channel_size = size(ChanMod_output{bb,uu}.H,3);
                    case 'FastFading'
                        channel_size = size(ChanMod_output{bb,uu}.H,4);
                end
                y_tx_filtered = zeros(size(BS_output(bb).y_tx,1) + channel_size - 1,ChanMod.nRX);
                ChanMod_output{bb,uu}.genie.channel_samples_position = nan(LTE_params.Nsub,ChanMod.nRX,ChanMod.nTX);
                for rr = 1:ChanMod.nRX
                    for tt = 1:ChanMod.nTX
                        switch ChanMod.filtering
                            case 'BlockFading'
                                y_tx_filtered(:,rr) = y_tx_filtered(:,rr) + conv(squeeze(ChanMod_output{bb,uu}.H(rr,tt,:)), BS_output(bb).y_tx(:,tt));
                            case 'FastFading'
                                taps_num = size(squeeze(ChanMod_output{bb,uu}.H(rr,tt,:,:)),2);    %number of channel taps
                                y_tx_filtered_matrix = [repmat(BS_output(bb).y_tx(:,tt),1,taps_num).*squeeze(ChanMod_output{bb,uu}.H(rr,tt,:,:))].'; %multiplication of the samples with their channel
                                %                         [y_tx_filtered_part mean_channel] = LTE_time_variant_convolution(y_tx_filtered_matrix.', squeeze(ChanMod_output.H(rr,tt,:,:)), LTE_params.Nsub, LTE_params.NfftCP{1}, LTE_params.NfftCP{2} , taps_num);
                                %                         spec = fft(mean_channel,LTE_params.Nfft,2);
                                if(strcmp(LTE_params.CyclicPrefix,'normal'))
                                    mapping_longer = logical(toeplitz([1;zeros(LTE_params.NfftCP{1}-1,1)],[ones(taps_num,1);zeros(LTE_params.NfftCP{1}-1,1)])); %the structur for time variant convolution for first and 7th symbols
                                    mapping_shorter = mapping_longer(1:LTE_params.NfftCP{2},1:(LTE_params.NfftCP{2} + taps_num - 1));   %mapping for other symbols, same as previous, but we just removed last column and row
                                    y_tx_filtered_part = zeros(length(BS_output(bb).y_tx)+taps_num-1,1);
                                    start = 1;
                                    % time-variant convolution
                                    %loop over ofdm symbols
                                    for symbol_i = 1:LTE_params.Nsub
                                        if(symbol_i == 1 || symbol_i == 7) %the calculation of the output for longer symbols
                                            y_tx_filtered_part_shifted_longer = zeros(LTE_params.NfftCP{1} + taps_num - 1,LTE_params.NfftCP{1});    %reseting of the convolution structure
                                            stop = start + LTE_params.NfftCP{1} + taps_num - 1 - 1;
                                            y_tx_filtered_part_shifted_longer(mapping_longer') = reshape(y_tx_filtered_matrix(:,start:stop - taps_num + 1),[],1);
                                            y_tx_filtered_part(start:stop,1) = y_tx_filtered_part(start:stop,1) + sum(y_tx_filtered_part_shifted_longer,2);
                                            %perfect channel is the mean of all channel realizations over one ofdm symbol
                                            %                                     spec = fft([mean(squeeze(ChanMod_output.H(rr,tt,start:stop - taps_num + 1,:)),1), zeros(1, LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,1,:))))]);
                                            %                                     ChanMod_output.genie.H_fft(:,symbol_i,rr,tt) = spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).'; % remove DC carrier and zeros padded up to size of FFT
                                            
                                            channel_time = zeros(LTE_params.NfftCP{1} + taps_num - 1,LTE_params.NfftCP{1});
                                            channel_time(mapping_longer') = reshape(squeeze(ChanMod_output{bb,uu}.H(rr,tt,start:stop - taps_num + 1,:)).',[],1);
                                            
                                            %mean channel position
                                            mean_channel = mean(squeeze(ChanMod_output{bb,uu}.H(rr,tt,start:stop - taps_num + 1,:)),1);
                                            [value mean_channel_position] = min(sum(abs(squeeze(ChanMod_output{bb,uu}.H(rr,tt,start:stop - taps_num + 1,:)) - repmat(mean_channel,LTE_params.NfftCP{1},1)),2));
                                            ChanMod_output{bb,uu}.genie.channel_samples_position(symbol_i,rr,tt) = mean_channel_position + start - 1;
                                            
                                            channel_time = channel_time(1:LTE_params.NfftCP{1},1:LTE_params.NfftCP{1});
                                            CP_length = LTE_params.NfftCP{1} - LTE_params.Nfft;
                                            F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                            F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                            DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                            channel_freq = DFT_matrix*F_CP_rem*channel_time*F_CP_add*DFT_matrix';
                                            ChanMod_output{bb,uu}.genie.H_fft_matrix(:,:,symbol_i,rr,tt) = channel_freq([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                            channel_freq_small = channel_freq([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                            ChanMod_output{bb,uu}.genie.H_fft(:,symbol_i,rr,tt) = channel_freq_small(diag(true(LTE_params.Ntot,1)));
                                            
                                            start = start + LTE_params.NfftCP{1};
                                        else %the calculation of the output for shorter symbols
                                            y_tx_filtered_part_shifted_shorter = zeros(LTE_params.NfftCP{2} + taps_num - 1,LTE_params.NfftCP{2});    %reseting of the convolution structure
                                            stop = start + LTE_params.NfftCP{2} + taps_num - 1 - 1;
                                            y_tx_filtered_part_shifted_shorter(mapping_shorter') = reshape(y_tx_filtered_matrix(:,start:stop - taps_num + 1),[],1);
                                            y_tx_filtered_part(start:stop,1) = y_tx_filtered_part(start:stop,1) + sum(y_tx_filtered_part_shifted_shorter,2);
                                            %perfect channel is the mean over all channel realizations over one ofdm symbol
                                            %                                     spec = fft([mean(squeeze(ChanMod_output.H(rr,tt,start:stop - taps_num + 1,:)),1), zeros(1, LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,1,:))))]);
                                            %                                     ChanMod_output.genie.H_fft(:,symbol_i,rr,tt) = spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).'; % remove DC carrier and zeros padded up to size of FFT
                                            
                                            channel_time = zeros(LTE_params.NfftCP{2} + taps_num - 1,LTE_params.NfftCP{2});
                                            channel_time(mapping_shorter') = reshape(squeeze(ChanMod_output{bb,uu}.H(rr,tt,start:stop - taps_num + 1,:)).',[],1);
                                            
                                            %mean channel position
                                            mean_channel = mean(squeeze(ChanMod_output{bb,uu}.H(rr,tt,start:stop - taps_num + 1,:)),1);
                                            [value mean_channel_position] = min(sum(abs(squeeze(ChanMod_output{bb,uu}.H(rr,tt,start:stop - taps_num + 1,:)) - repmat(mean_channel,LTE_params.NfftCP{2},1)),2));
                                            ChanMod_output{bb,uu}.genie.channel_samples_position(symbol_i,rr,tt) = mean_channel_position + start - 1;
                                            
                                            channel_time = channel_time(1:LTE_params.NfftCP{2},1:LTE_params.NfftCP{2});
                                            CP_length = LTE_params.NfftCP{2} - LTE_params.Nfft;
                                            F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                            F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                            DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                            channel_freq = DFT_matrix*F_CP_rem*channel_time*F_CP_add*DFT_matrix';
                                            ChanMod_output{bb,uu}.genie.H_fft_matrix(:,:,symbol_i,rr,tt) = channel_freq([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                            channel_freq_small = channel_freq([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                            ChanMod_output{bb,uu}.genie.H_fft(:,symbol_i,rr,tt) = channel_freq_small(diag(true(LTE_params.Ntot,1)));
                                            
                                            
                                            
                                            start = start + LTE_params.NfftCP{2};
                                        end
                                        %                             ChanMod_output.genie.H_fft(:,symbol_i,rr,tt) = spec(:,[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).';
                                    end
                                else
                                    mapping = logical(toeplitz([1;zeros(LTE_params.NfftCP-1,1)],[ones(taps_num,1);zeros(LTE_params.NfftCP-1,1)]));
                                    y_tx_filtered_part = zeros(length(BS_output.y_tx)+taps_num-1,1);
                                    start = 1;
                                    % time-variant convolution
                                    %loop over ofdm symbols
                                    for symbol_i = 1:LTE_params.Nsub
                                        y_tx_filtered_part_shifted_longer = zeros(LTE_params.NfftCP + taps_num - 1,LTE_params.NfftCP);    %reseting of the convolution structure
                                        stop = start + LTE_params.NfftCP + taps_num - 1 - 1;
                                        y_tx_filtered_part_shifted_longer(mapping') = reshape(y_tx_filtered_matrix(:,start:stop - taps_num + 1),[],1);
                                        y_tx_filtered_part(start:stop,1) = y_tx_filtered_part(start:stop,1) + sum(y_tx_filtered_part_shifted_longer,2);
                                        %perfect channel is the mean of all channel realizations over one ofdm symbol
                                        %                                 spec = fft([mean(squeeze(ChanMod_output.H(rr,tt,start:stop - taps_num + 1,:)),1), zeros(1, LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,1,:))))]);
                                        %                                 ChanMod_output.genie.H_fft(:,symbol_i,rr,tt) = spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).'; % remove DC carrier and zeros padded up to size of FFT
                                        start = start + LTE_params.NfftCP;
                                        %                                 ChanMod_output.genie.H_fft(:,symbol_i,rr,tt) = spec(:,[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).';
                                    end
                                end
                                %                         ChanMod_output.genie.H_fft(:,:,rr,tt) = spec(:,[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).'; % remove DC carrier and zeros padded up to size of FFT
                                y_tx_filtered(:,rr) = y_tx_filtered(:,rr) + y_tx_filtered_part;
                        end
                    end
                end
        end
        
        %% Introduce symbol timing offset, asynchronous multi eNodeB is possible
        if LTE_params.introduce_timing_offset
            if LTE_params.connection_table(bb,uu) % desired signals
                if (mod(subframe_i, 10)==1)&&(LTE_params.UE_config.timing_offset>0)
                    ChanMod_output{bb,uu}.y_rx = [zeros(LTE_params.UE_config.timing_offset, LTE_params.UE_config.nRX); y_tx_filtered];
                else
                    ChanMod_output{bb,uu}.y_rx = circshift(y_tx_filtered(1:length(BS_output(bb).y_tx),:), LTE_params.UE_config.timing_offset);
                end
            else % intefering signal
                if mod(subframe_i, 10)==1
                    ChanMod_output{bb,uu}.y_rx = [zeros(LTE_params.UE_config.timing_offset+LTE_params.Delay_diff, LTE_params.UE_config.nRX); y_tx_filtered];
                else
                    ChanMod_output{bb,uu}.y_rx = circshift(y_tx_filtered(1:length(BS_output(bb).y_tx),:), LTE_params.UE_config.timing_offset+LTE_params.Delay_diff);
                end
            end
        else
            ChanMod_output{bb,uu}.y_rx = y_tx_filtered(1:length(BS_output(bb).y_tx),:);
        end
        
        %% Mutlibase station addition        
        if LTE_params.connection_table(bb,uu)
            channel_gain = 1;
        else
            channel_gain = 10^(-LTE_params.Power_diff(uu)/10);
        end
        if bb == 1
            UE_input{uu}.input = channel_gain*ChanMod_output{bb,uu}.y_rx;
        else
            UE_input{uu}.input = UE_input{uu}.input  + channel_gain*ChanMod_output{bb,uu}.y_rx;
        end
       
    end
    
    %% Add zeros for SNR estimation
    if LTE_params.SNR_estimation 
        UE_input{uu}.input = [zeros(LTE_params.number_of_zeros_for_SNR_estimation,size(UE_input{uu}.input,2)); UE_input{uu}.input;zeros(LTE_params.number_of_zeros_for_SNR_estimation,size(UE_input{uu}.input,2))];
    end
    
    
    %% Add Noise
    noise_RandStream = LTE_params.noise_RandStream;
    n = 10^(-SNR/20)*(randn(noise_RandStream,size(UE_input{uu}.input)) + 1i*randn(noise_RandStream,size(UE_input{uu}.input)))/sqrt(2); %;/sqrt(ChanMod.nRX);
    v = sqrt(LTE_params.Nfft/LTE_params.Ntot) * n;
    UE_input{uu}.input = UE_input{uu}.input + v;
    
    %% Add genie information
    BS_output(bb).cell_genie.v{uu} = v;     %noise at RX-antennas (preFFT)
    BS_output(bb).cell_genie.n{uu} = n;     %noise after FFT at receiver (postFFT)
    
end
