function [channel_estimate] = LTE_channel_estimator_mini(LTE_params, rx_ref_symbols, RefSym, RefMapping, sigma_n2, y_rx_assembled)
% LTE channel estimator for mini receiver- to filter the output of the transmitter.
% [chan_output] = LTE_channel_model(BS_output, SNR)
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) by INTHFT
% www.nt.tuwien.ac.at
%
% input :   ChanMod                     ... struct -> include info about filtering type, channel autocorrelation matrix

% output:   channel_estimate            ... [ x ] matrix of estimate of channel coefficients
%
% date of creation: 2010/11/16

% global LTE_params;

%% Channel estimator
channel_estimate = nan(LTE_params.Ntot,LTE_params.Nsub,LTE_params.UE_config.nRX,LTE_params.BS_config.nTx);
interpolation = false;
channel_estimation_method = LTE_params.UE_config.channel_estimation_method;
channel_interpolation_method = LTE_params.UE_config.channel_interpolation_method;
filtering = 'BlockFading';

for tt = 1:LTE_params.BS_config.nTx
    for rr = 1:LTE_params.UE_config.nRX
        switch filtering
            case 'BlockFading'
                switch channel_estimation_method
                    case 'LS' %LS for block fading
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq),2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq)/2,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_one_channel = H_ls_one_channel(reference_mapping);
                        interpolation = true;

                    case  {'GENIE','GENIE_ALMMSE'}
                        %do nothing
                    otherwise
                        error('channel estimation not supported');
                end
                
                if(interpolation)
                    switch channel_interpolation_method
                        case {'linear', 'cubic', 'spline'}
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            channel_estimate_help = interp1(first_ref_sym:3:LTE_params.Ntot,H_ls_one_channel,1:LTE_params.Ntot,channel_interpolation_method,'extrap');
                            channel_estimate_help = channel_estimate_help.';
                        case 'sinc_freq'
                            pulse = zeros(LTE_params.Ntot,1);
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            pulse(first_ref_sym:3:end) = H_ls_one_channel;
                            channel_estimate_help = conv(pulse,sinc((-LTE_params.Ntot:LTE_params.Ntot)/3));
                            channel_estimate_help = channel_estimate_help((LTE_params.Ntot +1 ):(LTE_params.Ntot + LTE_params.Ntot)).';
                        case 'sinc_time'
                            pulse = zeros(LTE_params.Ntot,1);
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            pulse(first_ref_sym:3:end) = H_ls_one_channel;
                            PULSE = ifft(pulse,LTE_params.Nfft);
                            N = 3*[ones(floor(LTE_params.Nfft/(3*2)),1);zeros(LTE_params.Nfft - 2*floor(LTE_params.Nfft/(3*2)),1);ones(floor(LTE_params.Nfft/(3*2)),1)]; %rect
                            channel_time = PULSE.*N;
                            channel_estimate_help = fft(channel_time,LTE_params.Nfft);
                            channel_estimate_help = channel_estimate_help(1:LTE_params.Ntot);
                        case 'T-F'
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            channel_estimate_help = interpft(H_ls_one_channel,LTE_params.Ntot);
                            channel_estimate_help = circshift(channel_estimate_help,first_ref_sym-1);
                            
                        otherwise
                            error('interpolation method not supported')
                    end
                    channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help,1 , LTE_params.Nsub);
                end
            
            otherwise
                error('not supported filtering type')
        end
    end
end

switch filtering
    case 'BlockFading'
        switch channel_estimation_method
            case 'GENIE' %LS for block fading
                %                         A = [];
                %                         for ii = 1:LTE_params.Nsub
                %                             A = [A;eye(ChanMod.nRX)];
                %                         end
                %                         B = [];
                %                         for ii = 1:LTE_params.Nsub
                %                             B = [B,eye(ChanMod.nTX)];
                %                         end
                A = [];
                B = [];
                index_vec = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
                for sub_i = 1:LTE_params.Ntot
                    y_tx_assembled_genie_sub = conj(squeeze(y_tx_assembled_genie(sub_i,index_vec,:)));
                    if ChanMod.nTX == 1
                        y_tx_assembled_genie_sub = y_tx_assembled_genie_sub(:);
                    end
                    if rank(y_tx_assembled_genie_sub)~=ChanMod.nTX
                        A = [A;sub_i];
                    elseif rank(y_tx_assembled_genie_sub)==ChanMod.nTX
                        B = [B;sub_i];
                        y_rx_assembled_sub = conj(squeeze(y_rx_assembled(sub_i,index_vec,:)));
                        y_rx_assembled_sub = y_rx_assembled_sub(:);
                        kron_matrix = kron(eye(ChanMod.nRX),y_tx_assembled_genie_sub);
                        LS_matrix = inv(kron_matrix'*kron_matrix)*kron_matrix';
                        channel_help = LS_matrix*y_rx_assembled_sub;
                        channel_help = reshape(channel_help,ChanMod.nTX,ChanMod.nRX)';
                        channel_estimate(sub_i,1,:,:) = channel_help;
                    end
                    %y_tx_assembled_genie_sub =
                    %y_tx_assembled_genie_sub(:);
                end
                
                for tt = 1:ChanMod.nTX
                    for rr = 1:ChanMod.nRX
                        channel_estimate_help = interp1(B,squeeze(channel_estimate(B,1,rr,tt)),1:LTE_params.Ntot,channel_interpolation_method,'extrap');
                        channel_estimate_help = channel_estimate_help.';
                        channel_estimate(:,1,rr,tt) = channel_estimate_help;
                    end
                end
                channel_estimate = repmat(channel_estimate(:,1,:,:),[1,LTE_params.Nsub,1,1]);
        end
end