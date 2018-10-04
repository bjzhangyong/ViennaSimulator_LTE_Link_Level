function LTE_RX(LTE_params, ChanMod_output, UE_input, ChanMod, SNR, AtPort, subframe_i, BS, UE, UE_output, BS_output, uu)
% LTE receiver for a specific user.
% [UE_output, ACK, UE] = LTE_RX(chan_output, SNR, AtPort, subframe_i, BS, UE, BS_output)
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   chan_output         ... [1 x 1]struct - channel output:
%                                   [LTE_params.TxSymbols x 1]double receive signal y_rx
%           SNR                 ... [1 x 1]double current SNR value
%           AtPort              ... [1 x 1]double antenna port number
%           subframe_i          ... [1 x 1]double number of the subframe transmitted
%           BS                  ... [1 x nBS]struct - Base stations parameters from LTE_load_parameters.m
%           UE                  ... [1 x nUE]struct - User equipments capabilities from LTE_load_parameters.m
%           BS_output           ... [1 x nBS]struct - Base Station output
%                                   [LTE_params.TxSymbols x 1]double transmit signal y_tx,
%                                   [1 x nUE]struct - scheduler output and coding parameters, e.g:
%                                   [1 x # data bits]logical genie.data_bits
%                                   [1 x # sent bits(coded)]logical
%                                   genie.sent_bits
%                                   [1 x nUE]struct - user specific parameters and variables
% output:   UE_output           ... [1 x nUE]struct UE output for rx_coded_bits bits before decoder and after demapper
%                                                             and logical decoded bits
%           ACK                 ... [N_subframes x length(SNR_vec)]logical positive or negative acknowledgement
%           UE                  ... [1 x 1]struct - Updated UE struct
%
% date of creation: 2008/08/11
% last changes: 2008/09/08  Bosanska    Rewritten for the multi-user scenario
%               2008/09/10  Bosanska    SNR_i - new input variable for the correct track of ACK
%                                       ACK is now part of the UE_output structure
%               2008/09/15  Bosanska    changed structure of BS and UE (acc. to LTE_load_parameters.m)
%               2008/09/18  Bosanska    changed the name of chan_output -> ChanMod_output
%               2008/10/02  Bosanska    changed structure to multiple call of function LTE_RX
%                                       for multi-user scenario (parallel receivers)
%                                       added input [1x1]struct BS_UE_specific from BS_output.UE_specific(uu)
%                                       added input [1x1]double uu
%               2008/10/31  Bosanska    added matched receiver and perfect
%                                       channel knowledge for each user
%               2008/12/11  Simko       added disassambling of reference symbols
%                                       added channel estimator(LTE_channel_estimator)
% global Hsave

UE_signaling = BS_output.UE_signaling(uu);
genie_data = BS_output.genie(uu);

sigma_n2 = 10^(-SNR/10);
% Perform reception only if there is data to receive
if(BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs)
    
    % Get necessary data and convert in in the correct format
    nLayers    = BS_output.UE_signaling(uu).MCS_and_scheduling.nLayers;
    nCodewords = BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords;
    tx_mode    = BS_output.UE_signaling(uu).MCS_and_scheduling.tx_mode;
    
    % Code needed for legacy reasons
    %     switch BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords
    %         case 1
    %             UE_mapping = repmat(BS_output.UE_signaling(uu).MCS_and_scheduling.UE_mapping(:,1),1,2);
    %         case 2
    %             UE_mapping = BS_output.UE_signaling(uu).MCS_and_scheduling.UE_mapping;
    %         otherwise
    %             error('Number of codewords not defined');
    %     end
    UE_mapping = BS_output.UE_signaling(uu).MCS_and_scheduling.UE_mapping;
    %% Introduce carrier frequency offset
    if LTE_params.introduce_frequency_offset
        y_rx_unsync = UE_input.input.*exp(1i*2*pi*UE(uu).carrier_freq_offset*repmat((0:1:(size(ChanMod_output.y_rx,1)-1))',1,UE(uu).nRX)/LTE_params.Nfft);
    else
        y_rx_unsync = UE_input.input;
    end
    
    %% Estimate SNR
    if LTE_params.SNR_estimation 
        [S_plus_noise_power N_power] = LTE_SNR_estimator(LTE_params,y_rx_unsync);
        Noise_power = sum(mean(N_power,1));
        Signal_power = sum(S_plus_noise_power) - Noise_power;
        
        SNR = 10*log10(Signal_power/Noise_power);
        sigma_n2 = 10^(-SNR/10);
        y_rx_unsync(1:LTE_params.number_of_zeros_for_SNR_estimation,:) = [];
        %y_rx_unsync(end-LTE_params.number_of_zeros_for_SNR_estimation+1:end,:) = [];
        y_rx_unsync = y_rx_unsync(1:LTE_params.Fs*LTE_params.Tsubframe,:);
        
        UE_output(uu).Signal_plus_noise_power = S_plus_noise_power;
        UE_output(uu).Noise_power = N_power;
    else
        UE_output(uu).Signal_plus_noise_power = 1;
        UE_output(uu).Noise_power = sigma_n2;
    end
    
    %% Calculate the correct number of subframe from 1 - 10
    subframe_corr = mod(subframe_i,10);
    if(subframe_corr == 0)
        subframe_corr = 10;
    end
    
    %% Read pregenerated reference signals
    RefSym          = LTE_params.Reference_Signal(1,subframe_corr).RefSym;
    RefMapping      = LTE_params.Reference_Signal(1,subframe_corr).RefMapping;
    total_no_refsym = LTE_params.Reference_Signal(1,subframe_corr).total_no_refsym;
    NoData_indices  = LTE_params.Reference_Signal(1,subframe_corr).NoData_indices;
    % [RefSym,RefMapping,total_no_refsym,NoData_indices]=LTE_Common_gen_Reference_Signal(LTE_params.Nrb, LTE_params.Nsub, BS.nTX, BS.NIDcell, subframe_corr);
    
    %% Generation of synchronization signals and Initialization of signals
    rb_numbers = find(UE_mapping == 1);
    PrimMapping = zeros(size(NoData_indices));
    SecMapping = zeros(size(NoData_indices));
    CHmapping = zeros(size(NoData_indices));
    CHusedElements = zeros(LTE_params.Nrb*2,1);
    if LTE_params.usePBCH
        CHmapping = CHmapping + LTE_params.PBCH(1,subframe_corr).Mapping;
        CHusedElements = CHusedElements + LTE_params.PBCH(1,subframe_corr).UsedElements;
        %         CHnsyms_subframe = CHnsyms_subframe + LTE_params.PBCH(1,subframe_corr).N_Elements;
    end
    if LTE_params.usePDCCH
        CHmapping = CHmapping + LTE_params.PDCCH(1,subframe_corr).Mapping;
        CHusedElements = CHusedElements + LTE_params.PDCCH(1,subframe_corr).UsedElements;
        %         CHnsyms_subframe = CHnsyms_subframe + LTE_params.PDCCH(1,subframe_corr).N_Elements;
    end
    if(subframe_corr == 1 || subframe_corr == 6)
        % Read pregenerated sync signals
        PrimSync         = LTE_params.Sync_Signal(1,subframe_corr).PrimSync;
        PrimMapping      = LTE_params.Sync_Signal(1,subframe_corr).PrimMapping;
        SecSync          = LTE_params.Sync_Signal(1,subframe_corr).SecSync;
        SecMapping       = LTE_params.Sync_Signal(1,subframe_corr).SecMapping;
        SyncUsedElements = LTE_params.Sync_Signal(1,subframe_corr).SyncUsedElements;
        % [PrimSync, PrimMapping, SecSync, SecMapping, SyncUsedElements] = LTE_Common_gen_Synchronization_Signal(BS.NIDcell, LTE_params.Nrb, LTE_params.Nsc, LTE_params.Nsub, subframe_corr);
        
        %     NoData_indices = NoData_indices|PrimMapping;
        %     NoData_indices = NoData_indices|SecMapping;
        
        ResMapping = LTE_params.Sync_Signal(1,subframe_corr).ResMapping;    % reserved REs, nothing is transmitted here!
        CHusedElements = CHusedElements + LTE_params.Sync_Signal(1,subframe_corr).ResUsedElements;
        %         CHnsyms_subframe = CHnsyms_subframe + LTE_params.Sync_Signal(1,subframe_corr).ResNElements;
        CHmapping = CHmapping | ResMapping;
        
        rb_rx_symbols = LTE_params.Nsc*LTE_params.Ns - total_no_refsym - SyncUsedElements(UE_mapping) - CHusedElements(UE_mapping);
        rx_user_symbols = zeros(sum(rb_rx_symbols),UE(uu).nRX);
        H_est_user = zeros(sum(rb_rx_symbols),UE(uu).nRX,BS.nTX);
    elseif ~LTE_params.usePDCCH
        rx_symbols  = zeros((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
        PE_SINR_matrix = nan((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
        H_temp_user = zeros((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2,UE(uu).nRX,BS.nTX);
    else
        rb_rx_symbols = LTE_params.Nsc*LTE_params.Ns - total_no_refsym -  CHusedElements(UE_mapping);
        rx_user_symbols = zeros(sum(rb_rx_symbols),UE(uu).nRX);
        H_est_user = zeros(sum(rb_rx_symbols),UE(uu).nRX,BS.nTX);
    end
    
    %% Carrier Frequency Offset Compensation
    if (LTE_params.introduce_frequency_offset || LTE_params.introduce_timing_offset)
        if (subframe_corr == 1 || subframe_corr == 6)
            [y_rx_sync, UE_output(uu).freq_offset_est, UE_output(uu).timing_offset_est] = LTE_rx_sync(LTE_params, ChanMod, y_rx_unsync, UE(uu), UE_output(uu).freq_offset_est, UE_output(uu).timing_offset_est, subframe_corr, sigma_n2, RefSym, RefMapping, PrimSync, PrimMapping, SecSync, SecMapping);
        else
            [y_rx_sync, UE_output(uu).freq_offset_est, UE_output(uu).timing_offset_est] = LTE_rx_sync(LTE_params, ChanMod, y_rx_unsync, UE(uu), UE_output(uu).freq_offset_est, UE_output(uu).timing_offset_est, subframe_corr, sigma_n2, RefSym, RefMapping);
        end
        % calculate effective genie channel state information (for perfect channel knowledge)
        if ((UE_output(uu).timing_offset_est ~= UE.timing_offset) && LTE_params.introduce_timing_offset)%&&(strcmp(UE.channel_estimation_method,'PERFECT')))
            for rr = 1:ChanMod.nRX
                for tt = 1:ChanMod.nTX
                    spec = fft(circshift([squeeze(ChanMod_output.H(rr,tt,:)); zeros(LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,:))),1)],UE.timing_offset-UE_output(uu).timing_offset_est));
                    ChanMod_output.genie.H_fft(:,:,rr,tt) = repmat(spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]),1,LTE_params.Nsub); % remove DC carrier and zeros padded up to size of FFT
%                     figure
%                     plot(abs(spec))
%                     figure
%                     plot(angle(spec))
                end
            end
        end
    else
        y_rx_sync = y_rx_unsync;
        UE_output(uu).freq_offset_est.error = NaN;
        UE_output(uu).freq_offset_est.frac = NaN;
        UE_output(uu).freq_offset_est.int = NaN;
        UE_output(uu).freq_offset_est.res = NaN;
    end
    
    
    %% Remove CP, FFT, remove zeros
    for nn = 1:UE(uu).nRX
        y_rx_resolved = cell(4,1);
        if(length(LTE_params.Tg)==2)
            y_rx_resolved{1} = y_rx_sync(1:LTE_params.NfftCP{1},nn);
            y_rx_resolved{2} = reshape(y_rx_sync(LTE_params.NfftCP{1}+1:LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
            y_rx_resolved{3} = y_rx_sync(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn);
            y_rx_resolved{4} = reshape(y_rx_sync(2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:end,nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
            y_rx_assembled_ifft = [y_rx_resolved{1}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{2}(LTE_params.Index_RxCyclicPrefix{2},:)...
                y_rx_resolved{3}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{4}(LTE_params.Index_RxCyclicPrefix{2},:)];
        else
            y_rx = reshape(y_rx_sync(:,nn),LTE_params.NfftCP,LTE_params.Nsub);
            y_rx_assembled_ifft = y_rx(LTE_params.Index_RxCyclicPrefix,:);
        end
        
        y_rx_assembled_shifted = fft(1/sqrt(LTE_params.Nfft)/(sqrt(LTE_params.Nfft/LTE_params.Ntot))*y_rx_assembled_ifft);
        y_rx_assembled_padded = circshift(y_rx_assembled_shifted,LTE_params.Ntot/2);
        % remove zero DC carrier
        y_rx_assembled(:,:,nn) = y_rx_assembled_padded([1:LTE_params.Ntot/2,LTE_params.Ntot/2+2:LTE_params.Ntot+1],:);
    end
    %     power = [power,mean(mean(abs(y_rx_assembled).^2,1),2)];
    %% Disassemble reference symbols
    rx_ref_symbols = nan([size(RefSym),ChanMod.nRX]);   %allocate memery for reference signal for every channel, the number of transmitt antennas is included in RefSym
    for tt = 1:ChanMod.nTX
        for rr = 1:ChanMod.nRX
            rx_ref_symbols_help = y_rx_assembled(:,:,rr);   %use recieved symbols from one recieve antenna
            rx_ref_symbols_help = rx_ref_symbols_help(RefMapping(:,:,tt));  %extract the signal on pilots positons
            if(tt>2)    %in case of 3rd and 4th transmitt antenna, there are less pilots, so we have to insert some nan symbols,to be consistent
                rx_ref_symbols_help = [rx_ref_symbols_help;nan(size(rx_ref_symbols_help))];
            end
            rx_ref_symbols(:,:,tt,rr) = reshape(rx_ref_symbols_help,size(RefSym(:,:,tt)));  %finally place the reference signal on allocated position
        end
    end
    
    %% Channel and noise estimation
    H_est = LTE_channel_estimator(LTE_params,ChanMod, UE(uu), rx_ref_symbols, RefSym, RefMapping, ChanMod_output.genie.H_fft,sigma_n2,BS_output.cell_genie.y_tx_assembled,y_rx_assembled);
    %     H_est = ChanMod_output.genie.H_fft;
    if isempty(ChanMod_output.genie.H_fft)
        ChanMod_output.genie.H_fft = zeros(LTE_params.Ntot,LTE_params.Nsub,LTE_params.UE_config.nRX,LTE_params.BS_config.nTx);
    end
    if ~strcmp('TB',LTE_params.Simulation_type)
        UE_output(uu).channel_estimation_error = reshape(sum(sum(abs(ChanMod_output.genie.H_fft - H_est).^2,1),2)/(LTE_params.Nsub*LTE_params.Ntot),size(H_est,3),size(H_est,4));
    else
        UE_output(uu).channel_estimation_error = 0;
    end
    %     channel_estimation_error_freq_depend(:,:,:,:,counti) =  abs(ChanMod_output.genie.H_fft(:,:,:,:) - H_est(:,:,:,:)).^2;
    %     counti = counti+1;
    if strcmp('ICI_aware_ZF',UE(uu).receiver)
        ChanMod_output.genie.H_fft_matrix_est = LTE_ICI_estimator(LTE_params,H_est,ChanMod,ChanMod_output);
    end
    %% Disassemble symbols
    
    for nn = 1:UE(uu).nRX
        if(subframe_corr == 1 || subframe_corr == 6 || LTE_params.usePDCCH)
            index = 0;
            for ii = 1:length(rb_numbers)
                if(rb_numbers(ii) > LTE_params.Nrb)
                    y_rx_assembled_temp = y_rx_assembled((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end,nn);
                    H_temp =                reshape(H_est((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end,nn,:),[LTE_params.Nsc*LTE_params.Ns,1,LTE_params.BS_config.nTx]);
                    rx_user_symbols(index+1:index+rb_rx_symbols(ii),nn) = y_rx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1-LTE_params.Nrb)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,8:8+LTE_params.Ns-1)));
                    H_est_user(index+1:index+rb_rx_symbols(ii),nn,:) =                 H_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1-LTE_params.Nrb)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,8:8+LTE_params.Ns-1)),:,:);
                else
                    y_rx_assembled_temp = y_rx_assembled((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns,nn);
                    H_temp =                reshape(H_est((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns,nn,:),[LTE_params.Nsc*LTE_params.Ns,1,LTE_params.BS_config.nTx]);
                    rx_user_symbols(index+1:index+rb_rx_symbols(ii),nn) = y_rx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)...
                        + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)));
                    H_est_user(index+1:index+rb_rx_symbols(ii),nn,:) =                 H_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)...
                        + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)),:,:);
                end
                index = index + rb_rx_symbols(ii);
            end
        else
            for ii = 1:LTE_params.Nrb
                for tt = 1:BS.nTX
                    H_temp =                H_est((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,1:LTE_params.Ns,nn,tt);
                    H_temp_user(:,ii,1,nn,tt) =      H_temp(~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                    H_temp =                H_est((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,LTE_params.Ns+1:end,nn,tt);
                    H_temp_user(:,ii,2,nn,tt) =      H_temp(~NoData_indices(1:LTE_params.Nsc,LTE_params.Ns + 1:end));
                end
                % NOTE: some more comments
                y_rx_assembled_temp = y_rx_assembled((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,1:LTE_params.Ns,nn);
                rx_symbols(:,ii,1) = y_rx_assembled_temp(~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                y_rx_assembled_temp = y_rx_assembled((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,LTE_params.Ns+1:end,nn);
                rx_symbols(:,ii,2) = y_rx_assembled_temp(~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                %             H_temp_user(ii,2,:,nn) = H_temp_user(ii,1,:,nn);
            end
            %rx_symbols_temp = rx_symbols(:,:,:);
            rx_user_symbols(:,nn) = rx_symbols(repmat(reshape(UE_mapping,[1 size(UE_mapping)]),[LTE_params.Nsc*LTE_params.Ns - ...
                sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))),1,1]));
            for tt = 1:BS.nTX
                H_temp_user_temp = squeeze(H_temp_user(:,:,:,nn,tt));
                H_est_user(:,nn,tt) =   H_temp_user_temp(repmat(reshape(UE_mapping,[1 size(UE_mapping)]),[LTE_params.Nsc*LTE_params.Ns - ...
                    sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))),1,1]));
            end
        end
    end
    
    %% Perform detection
    %     MSE = LTE_channelestimator_MSE(sigma_n2,UE(uu).channel_autocorrelation_matrix,LTE_params,BS.nAtPort,UE(uu).user_speed);
    MSE = zeros(LTE_params.Ntot,1);
    switch LTE_params.UE_config.receiver
        case 'ICI_aware_ZF'
            if(subframe_corr == 1 || subframe_corr == 6)
                [LLR_SD,M,rx_layer_x_equalized interference_power signal_power] = LTE_ici_aware_detecting(LTE_params,UE_signaling,rx_user_symbols,subframe_corr,ChanMod,rb_numbers,y_rx_assembled,BS_output.genie(uu).layer_x,ChanMod_output.genie.H_fft_matrix_est,uu,RefSym,RefMapping,NoData_indices,PrimMapping,SecMapping,PrimSync,SecSync,CHmapping,rb_rx_symbols);
            else
                [LLR_SD,M,rx_layer_x_equalized interference_power signal_power] = LTE_ici_aware_detecting(LTE_params,UE_signaling,rx_user_symbols,subframe_corr,ChanMod,rb_numbers,y_rx_assembled,BS_output.genie(uu).layer_x,ChanMod_output.genie.H_fft_matrix_est,uu,RefSym,RefMapping,NoData_indices);
            end
        otherwise
            [LLR_SD,M,H_back,rx_layer_x_equalized] = LTE_detecting(UE_signaling.MCS_and_scheduling,BS.nAtPort,rx_user_symbols,UE(uu).nRX,H_est_user,ChanMod.filtering,H_est,LTE_params,UE(uu).receiver,UE(uu).receiver_k,sigma_n2,MSE);
    end
    %% SINR Calculation
    
    PE_noise_power_overall = abs(rx_layer_x_equalized-BS_output.genie(uu).layer_x.').^2;
    PE_signal_power_overall = abs(rx_layer_x_equalized).^2;
    if UE(uu).mode == 2
        PE_noise_power_overall = PE_noise_power_overall.';
        PE_noise_power_overall = PE_noise_power_overall(:);
        
        PE_signal_power_overall = PE_signal_power_overall.';
        PE_signal_power_overall = PE_signal_power_overall(:);
    end
    PE_noise_power_overall = mean(PE_noise_power_overall,2);
    PE_signal_power_overall = mean(PE_signal_power_overall,2);
    
    if(subframe_corr == 1 || subframe_corr == 6 || LTE_params.usePDCCH)
        index = 0;
        for ii = 1:length(rb_numbers)
            
            if(rb_numbers(ii) > LTE_params.Nrb)
                
                PE_noise_power_temp_vec = PE_noise_power_overall(index+1:index+rb_rx_symbols(ii),:);
                PE_noise_power_temp = nan(LTE_params.Nsc,LTE_params.Ns);
                PE_noise_power_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1-LTE_params.Nrb)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,8:8+LTE_params.Ns-1))) = PE_noise_power_temp_vec;
                UE_output(uu).PE_noise_power_subframe((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end) = PE_noise_power_temp;
                
                PE_signal_power_temp_vec = PE_signal_power_overall(index+1:index+rb_rx_symbols(ii),:);
                PE_signal_power_temp = nan(LTE_params.Nsc,LTE_params.Ns);
                PE_signal_power_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1-LTE_params.Nrb)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,8:8+LTE_params.Ns-1))) = PE_signal_power_temp_vec;
                UE_output(uu).PE_signal_power_subframe((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end) = PE_signal_power_temp;
                
            else
                
                PE_noise_power_temp_vec = PE_noise_power_overall(index+1:index+rb_rx_symbols(ii),:);
                PE_noise_power_temp = nan(LTE_params.Nsc,LTE_params.Ns);
                PE_noise_power_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns))) = PE_noise_power_temp_vec;
                UE_output(uu).PE_noise_power_subframe((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) = PE_noise_power_temp;
                
                PE_signal_power_temp_vec = PE_signal_power_overall(index+1:index+rb_rx_symbols(ii),:);
                PE_signal_power_temp = nan(LTE_params.Nsc,LTE_params.Ns);
                PE_signal_power_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns))) = PE_signal_power_temp_vec;
                UE_output(uu).PE_signal_power_subframe((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) = PE_signal_power_temp;
                
            end
            index = index + rb_rx_symbols(ii);
        end
    else
        UE_mapping_temp = UE_mapping(:);
        nr_of_data = LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)));
        index = 1;
        for ii = 1:LTE_params.Nrb*2
            if UE_mapping_temp(ii)
                PE_noise_power_temp_vec = PE_noise_power_overall(index:index+nr_of_data-1);
                PE_noise_power_temp = nan(LTE_params.Nsc,LTE_params.Ns);
                PE_noise_power_temp(~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)) = PE_noise_power_temp_vec;
                
                PE_signal_power_temp_vec = PE_signal_power_overall(index:index+nr_of_data-1);
                PE_signal_power_temp = nan(LTE_params.Nsc,LTE_params.Ns);
                PE_signal_power_temp(~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)) = PE_signal_power_temp_vec;
                
                if(ii > LTE_params.Nrb)
                    UE_output(uu).PE_noise_power_subframe((ii-LTE_params.Nrb-1)*LTE_params.Nsc+1:(ii-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end) = PE_noise_power_temp;
                    UE_output(uu).PE_signal_power_subframe((ii-LTE_params.Nrb-1)*LTE_params.Nsc+1:(ii-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end) = PE_signal_power_temp;
                else
                    UE_output(uu).PE_noise_power_subframe((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,1:LTE_params.Ns) = PE_noise_power_temp;
                    UE_output(uu).PE_signal_power_subframe((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,1:LTE_params.Ns) = PE_signal_power_temp;
                end
                index = index + nr_of_data;
            end
            
        end
    end
    
    %% Undo layer mapping
    % So we have codewords again
    switch BS.nAtPort
        case 1  % single antenna transmission
            LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
        case 2
            if (UE_signaling.MCS_and_scheduling.tx_mode == 2)     % transmit diversity
                LLR_SS{1} = reshape([LLR_SD(M(1):-1:1,:);LLR_SD(M(1)+M(2):-1:M(1)+1,:)],1,[]).';
            else           % spatial multiplexing
                switch nLayers
                    case 1
                        LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                    case 2
                        if (UE_signaling.MCS_and_scheduling.nCodewords == 2)
                            LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                            LLR_SS{2} = reshape(LLR_SD(M(1)+M(2):-1:M(1)+1,:),1,[]).';
                        else
                            LLR_SS{1} = reshape([LLR_SD(M(1):-1:1,:);LLR_SD(M(1)+M(2):-1:M(1)+1,:)],1,[]).';
                        end
                end
            end
        case 4
            if (UE_signaling.MCS_and_scheduling.tx_mode == 2)     % transmit diversity
                %             rx_user_symbols_mrc =
                %             reshape([rx_layer_x(1,:);conj(rx_layer_x(2,:));rx_layer_x(3,:);conj(rx_layer_x(4,:))],1,[]); %if ZF Detector is used then use this
                LLR_SS{1} = reshape([LLR_SD(M(1):-1:1,:);LLR_SD(M(1)+M(2):-1:M(1)+1,:);LLR_SD(M(1)+M(2)+M(3):-1:M(1)+M(2)+1,:);LLR_SD(M(1)+M(2)+M(3)+M(4):-1:M(1)+M(2)+M(3)+1,:)],1,[]).';
            else           % spatial multiplexing
                switch nLayers
                    case 1
                        LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                        %                     rx_symbols{1}=received_layer_x;
                    case 2
                        if (UE_signaling.MCS_and_scheduling.nCodewords == 2)
                            LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                            %                         LLR_SS{2} = reshape(LLR_SD(2*M(2):-1:M(2)+1,:),1,[]).';
                            LLR_SS{2} = reshape(LLR_SD(M(1)+M(2):-1:M(1)+1,:),1,[]).';
                        else
                            %                         LLR_SS{1} = reshape([LLR_SD(M(1):-1:1,:);LLR_SD(2*M(2):-1:M(2)+1,:)],1,[]).';
                            LLR_SS{1} = reshape([LLR_SD(M(1):-1:1,:);LLR_SD(M(1)+M(2):-1:M(1)+1,:)],1,[]).';
                        end
                    case 3
                        LLR_SS{1} = reshape(LLR_SD(M(1):-1:1,:),1,[]).';
                        %                     LLR_SS{2}= reshape(LLR_SD([2*M(2):-1:M(2)+1 3*M(2):-1:2*M(2)+1],:),1,[]).';
                        LLR_SS{2}= reshape(LLR_SD([M(1)+M(2):-1:M(1)+1 M(1)+M(2)+M(3):-1:M(1)+M(2)+1],:),1,[]).';
                    case 4
                        %                     LLR_SS{1} = reshape(LLR_SD([M(1):-1:+1 2*M(1):-1:M(1)+1],:),1,[]).';
                        %                     LLR_SS{2} = reshape(LLR_SD([3*M(2):-1:2*M(2)+1 4*M(2):-1:3*M(2)+1],:),1,[]).';
                        LLR_SS{1} = reshape(LLR_SD([M(1):-1:+1 M(1)+M(2):-1:M(1)+1],:),1,[]).';
                        LLR_SS{2} = reshape(LLR_SD([M(3)+M(2)+M(1):-1:M(2)+M(1)+1 M(4)+M(3)+M(2)+M(1):-1:M(3)+M(2)+M(1)+1],:),1,[]).';
                end
            end
    end
    
    %% Decoding
    for i = 1:UE_signaling.MCS_and_scheduling.nCodewords
        %%%%%%%%%%%%% Soft Sphere Decoder Solution
        rx_scrambled_bits{i} =  (1+sign(LLR_SS{i}.'))/2; % NOTE: I don't get what this line is for. Is to get the uncoded BER performance? => yup, exactly, and uncoded throughput and so on
        %% Decrambling of the bits, as of TS 36.211 V8.2.0 (2008-03) Section 6.3.1
        % for rx_coded_bits there is mode 'scramble', because descrambler operates on a bit basis that is equal to scrambling operation
        UE_output(uu).rx_coded_bits{i} = LTE_common_scrambling(rx_scrambled_bits{i},BS.NIDcell,uu,subframe_corr,i,'scramble');
        %UE_output(uu).rx_coded_bits{i} = rx_scrambled_bits{i};
        % for rx_data_bits the LLRs have to be descrambled first
        LLR_SS_descrambled{i} = LTE_common_scrambling(LLR_SS{i},BS.NIDcell,uu,subframe_corr,i,'descramble');
        %     LLR_SS_descrambled{i} = LLR_SS{i};
        %% DL-SCH: Decoding of the bits
        [UE_output(uu).rx_data_bits{i} BER_error_bits_turbo UE_output(uu).ACK(i) UE_output(uu).ACK_codeblocks(i) UE_output(uu).C(i)] = LTE_rx_DLSCH_decode(LTE_params,LLR_SS_descrambled{i},UE_signaling,UE(uu),BS.UE_specific(uu),genie_data,i);
        %         %% decrease the traffic model data buffer if successfully transmitted
        %         UE(uu).traffic_model.decrease_data_buffer(length(UE_output(uu).rx_data_bits{i}),UE_output(uu).ACK(i));
    end
    
    % Feedback performed only when data is received
    feedback_rv_idx = zeros(1,length(UE_signaling.turbo_rate_matcher));
    for rv_i = 1:length(UE_signaling.turbo_rate_matcher)
        feedback_rv_idx(rv_i) = UE_signaling.turbo_rate_matcher(rv_i).rv_idx;
    end
    UE_output(uu).rv_idx = feedback_rv_idx;
    % Precoding feedback calculation
    if LTE_params.uplink_delay ~= 0
        %             H_est = LTE_channel_estimator(LTE_params,ChanMod, UE(uu), rx_ref_symbols, RefSym, RefMapping, ChanMod_output.genie.H_fft,sigma_n2,BS_output.cell_genie.y_tx_assembled,y_rx_assembled);
        %             [RI_tmp,PMI_tmp]=LTE_feedback_precoding(BS.nAtPort,sigma_n2,LTE_params,UE_signaling.MCS_and_scheduling(1).CQI_params(1).modulation_order(1),ChanMod_output.genie.H_fft,UE(uu),uu);
        %             if LTE_params.UE_config.mode == 3
        %                 [RI_tmp,PMI_tmp,CQI_tmp]=LTE_feedback(BS.nAtPort,sigma_n2,LTE_params,H_back,UE(uu),uu,LTE_params.UE_config.mode);
        %             else
        [RI_tmp,PMI_tmp,CQI_tmp]=LTE_feedback(BS.nAtPort,sigma_n2,LTE_params,H_est,UE(uu),uu,LTE_params.UE_config.mode);
        %             end
        UE_output(uu).RI = RI_tmp;
        UE_output(uu).PMI = PMI_tmp-1;
        UE_output(uu).CQI = CQI_tmp;
    end
    
else
    % Feedback for the cases where no data was transmitted
    UE_output(uu).rv_idx = zeros(1,2);
    UE_output(uu).ACK(:) = 0;
    H_est = ChanMod_output.genie.H_fft;
    % Precoding feedback calculation. Use QPSK modulation for calculations when no data was assigned (assume the channel was really bad)
    if LTE_params.uplink_delay ~= 0
        [RI_tmp,PMI_tmp,CQI_tmp]=LTE_feedback(BS.nAtPort,sigma_n2,LTE_params,H_est,UE(uu),uu,LTE_params.UE_config.mode);
        UE_output(uu).RI = RI_tmp;
        UE_output(uu).PMI = PMI_tmp-1;
        UE_output(uu).CQI = CQI_tmp;
    end
    
end

% Finish feedback construction
% For now use perfect channel information. Can be changed to an SINR computed from channel estimation
% SINR = BS_output.cell_genie.SINR(uu,:);
% SINRs_to_CQI = SINR(1:12:end);
% UE_output(uu).CQI_feedback = floor(LTE_common_CQI_mapping(LTE_params.CQI_mapping,SINRs_to_CQI)); % only one CQI value over bandwidth, because there is only noise added
% UE_output(uu).CQI = UE_output(uu).CQI_feedback;
% UE_output(uu).PMI = 1;
% HARQ process id of this current ACK
UE_output(uu).HARQ_process = BS_output.UE_signaling(uu).MCS_and_scheduling.HARQ_process_id;
% [BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs]

UE_output(uu).UE_scheduled = repmat((BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs>0),BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords);

%% Calculate noise on a subcarrier basis with the genie information
% noise_SC = BS_output.cell_genie.y_tx_assembled - y_rx_assembled;
% SINR_SC_linear = mean(abs(noise_SC).^2,2);
% BS_output.cell_genie.genie_SNR = SINR_SC_linear;
