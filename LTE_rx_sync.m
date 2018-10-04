function [y_rx_sync, freq_offset_est, timing_offset_est] = LTE_rx_sync(LTE_params, ChanMod, y_rx, UE, freq_offset_est, timing_offset_est, subframe_corr, sigma_n2, RefSym, RefMapping, PrimSync, PrimMapping, SecSync, SecMapping)
% LTE symbol timing and carrier frequency synchronization
% Author: Qi Wang, qwang@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   y_rx                        ... [ XX x nRX ] received sequence
%           UE                          ... UE parameter
%           freq_offset_est             ... struct
%           RefSym                      ... Reference Symbol
%           RefMapping                  ... Reference Symbol Mapping
%           PrimSync                    ... Primary Synchronization signals
%           PrimMapping                 ... Primary Synchronization signals mapping
%           SecSync                     ... Secondary Synchronization signals
%           SecMapping                  ... Secondary Synchronization signals mapping
% output:   y_rx_sync                   ... [ XX x nRX ] sequence with carrier frequency offset compensated
%           freq_offset_est             ... struct
%
% date of creation: 2009/04/24
% last changes: 2009/06/23 Wang

if LTE_params.introduce_timing_offset
    switch UE.timing_sync_method
        case 'perfect'
            timing_offset_est = UE.timing_offset;
        case 'none'
            timing_offset_est = NaN;
        case 'autocorrelation'
            % J. J. van de Beek, M. Sandell, P. O. Borjesson, "ML
            % Estimation of Time and Frequency Offset in OFDM Systems",
            % IEEE Transactions on Signal processing, vol. 45, No. 7, July
            % 1997
            if (subframe_corr == 1)||(subframe_corr == 6)
                search_range = LTE_params.NfftCP{1}+LTE_params.Nfft;
                vv = zeros(search_range,1);
                if (length(LTE_params.Ng) == 2)
                    Index_extractCP_temp{1} = 1:LTE_params.Ng(1);
                    Index_extractCP_temp{2} = reshape(LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
                    Index_extractCP_temp{3} = LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+(1:LTE_params.Ng(1));
                    Index_extractCP_temp{4} = reshape(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-2)),LTE_params.NfftCP{2},(LTE_params.Ns-2));
                    Index_extractCP = [Index_extractCP_temp{1} reshape(Index_extractCP_temp{2}(1:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3} reshape(Index_extractCP_temp{4}(1:LTE_params.Ng(2),:),1,[])];
                    %             Index_extractCP = [Index_extractCP_temp{1}(3:end) reshape(Index_extractCP_temp{2}(3:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(3:end) reshape(Index_extractCP_temp{4}(3:LTE_params.Ng(2),:),1,[])];
                else
                    % undebugged
                    Index_extractCP_temp = reshape((1:LTE_params.NfftCP*LTE_params.Nsub),LTE_params.NfftCP,LTE_params.Nsub);
                    Index_extractCP = reshape(Index_extractCP_temp(1:LTE_params.Ng,:),1,[]);
                end
                for search_point = 1:search_range
                    vv(search_point) = sum(sum(y_rx(search_point-1+Index_extractCP,:).*conj(y_rx(search_point-1+Index_extractCP+LTE_params.Nfft,:))));
                end
                [~, ss] = max(vv);
                timing_offset_est = ss-1;
            end
        case 'min_energy'
            % M. Speth, F. Classen, H. Meyr, "Frame Synchronization of OFDM Systems in Frequency Selective Fading Channels", IEEE 47th VTC 1997. 
            if (subframe_corr == 1)||(subframe_corr == 6)
                search_range = LTE_params.NfftCP{2};
                vv = zeros(search_range,1);
                if (length(LTE_params.Ng) == 2)
                    Index_extractCP_temp{1} = 1:LTE_params.Ng(1);
                    Index_extractCP_temp{2} = reshape(LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
                    Index_extractCP_temp{3} = LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+(1:LTE_params.Ng(1));
                    Index_extractCP_temp{4} = reshape(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-2)),LTE_params.NfftCP{2},(LTE_params.Ns-2));
                    Index_extractCP = [Index_extractCP_temp{1}(end) reshape(Index_extractCP_temp{2}(LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(end) reshape(Index_extractCP_temp{4}(LTE_params.Ng(2),:),1,[])];
                    %             Index_extractCP = [Index_extractCP_temp{1}(3:end) reshape(Index_extractCP_temp{2}(3:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(3:end) reshape(Index_extractCP_temp{4}(3:LTE_params.Ng(2),:),1,[])];
                else
                    % undebugged
                    Index_extractCP_temp = reshape((1:LTE_params.NfftCP*LTE_params.Nsub),LTE_params.NfftCP,LTE_params.Nsub);
                    Index_extractCP = reshape(Index_extractCP_temp(1:LTE_params.Ng,:),1,[]);
                end
                for search_point = 1:search_range
                    metric_d = mean((abs(abs(y_rx(search_point-1+Index_extractCP+LTE_params.Nfft,:))-abs(y_rx(search_point-1+Index_extractCP,:)))).^2);
                    vv(search_point) = metric_d;
                end
                [~, ss] = min(vv);
                timing_offset_est = ss-1;
            end
        case 'energy_transition'
            % V. Le Nir, T. van Waterschoot,J. Duplicy, M. Moonen, "Blind coarse timing offset estimation for CP-OFDM and
            % ZP-OFDM transmission over frequency selective channels", EURASIP Journal on Wireless Communications and Networking
            % Volume 2009 (2009), Article ID 262813, 8 pages
            % doi:10.1155/2009/262813
            if (subframe_corr == 1)||(subframe_corr == 6)
                search_range = LTE_params.NfftCP{2};
                metric_d = zeros(LTE_params.Nsub-1,2);
                vv = zeros(search_range,1);
                if (length(LTE_params.Ng) == 2)
                    Index_extractCP_temp{1} = 1:LTE_params.Ng(1);
                    Index_extractCP_temp{2} = reshape(LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
                    Index_extractCP_temp{3} = LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+(1:LTE_params.Ng(1));
                    Index_extractCP_temp{4} = reshape(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-2)),LTE_params.NfftCP{2},(LTE_params.Ns-2));
                    Index_extractCP = [Index_extractCP_temp{1}(end) reshape(Index_extractCP_temp{2}(LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(end) reshape(Index_extractCP_temp{4}(LTE_params.Ng(2),:),1,[])];
                    %             Index_extractCP = [Index_extractCP_temp{1}(3:end) reshape(Index_extractCP_temp{2}(3:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(3:end) reshape(Index_extractCP_temp{4}(3:LTE_params.Ng(2),:),1,[])];
                else
                    % undebugged
                    Index_extractCP_temp = reshape((1:LTE_params.NfftCP*LTE_params.Nsub),LTE_params.NfftCP,LTE_params.Nsub);
                    Index_extractCP = reshape(Index_extractCP_temp(1:LTE_params.Ng,:),1,[]);
                end
                for search_point = 1:search_range
                    for kk = 0:1
                        metric_d(:,kk+1) = mean((abs(abs(y_rx(search_point-1+Index_extractCP+kk+LTE_params.Nfft,:))-abs(y_rx(search_point-1+Index_extractCP+kk,:)))).^2,2);
                    end
                    vv_temp = mean(metric_d);
                    vv(search_point) = sum(vv_temp(2))./sum(vv_temp(1));
                    
                end
                [~, ss] = max(vv);
                timing_offset_est = ss-1;
            end
        case 'modified_energy_transition'
            if (subframe_corr == 1)||(subframe_corr == 6)
                index_k = LTE_params.Ng(2)*2;
                search_range = LTE_params.NfftCP{2}-index_k/2;
                metric_d = zeros(LTE_params.Nsub-1,index_k);
                vv = zeros(search_range,1);
                if (length(LTE_params.Ng) == 2)
                    Index_extractCP_temp{1} = 1:LTE_params.Ng(1);
                    Index_extractCP_temp{2} = reshape(LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
                    Index_extractCP_temp{3} = LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+(1:LTE_params.Ng(1));
                    Index_extractCP_temp{4} = reshape(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-2)),LTE_params.NfftCP{2},(LTE_params.Ns-2));
                    Index_extractCP = [Index_extractCP_temp{1}(end) reshape(Index_extractCP_temp{2}(LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(end) reshape(Index_extractCP_temp{4}(LTE_params.Ng(2),:),1,[])];
                    %             Index_extractCP = [Index_extractCP_temp{1}(3:end) reshape(Index_extractCP_temp{2}(3:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(3:end) reshape(Index_extractCP_temp{4}(3:LTE_params.Ng(2),:),1,[])];
                else
                    % undebugged
                    Index_extractCP_temp = reshape((1:LTE_params.NfftCP*LTE_params.Nsub),LTE_params.NfftCP,LTE_params.Nsub);
                    Index_extractCP = reshape(Index_extractCP_temp(1:LTE_params.Ng,:),1,[]);
                end
                % estimate maximum channel delay
                for search_point = 1:search_range
                    for kk = -index_k/2+1:index_k/2
                        metric_d(:,kk+index_k/2) = sum((abs(abs(y_rx(search_point-1+Index_extractCP+kk+LTE_params.Nfft,:))-abs(y_rx(search_point-1+Index_extractCP+kk,:)))).^2,2);
                    end
                    vv_temp = sum(metric_d);
                    vv(search_point) = log(sum(vv_temp(index_k/2+1:end)))-log(sum(vv_temp(1:index_k/2)));
                end
                [~, ss] = sort(abs(vv));
                L = abs(ss(end)-ss(end-1));
                if L>LTE_params.Ng(1)
                    index_k = (LTE_params.Ng(1)-1)*2;
                else
                    index_k = (LTE_params.Ng(1)-L)*2;
                end
                % estimate the symbol timing
                search_range = LTE_params.NfftCP{2}-index_k/2;
                vv = zeros(search_range,1);
                metric_d = zeros(LTE_params.Nsub-1,index_k);
                for search_point = 1:search_range
                    for kk = -index_k/2+1:index_k/2
                        metric_d(:,kk+index_k/2) = mean((abs(abs(y_rx(search_point-1+Index_extractCP+kk+LTE_params.Nfft,:))-abs(y_rx(search_point-1+Index_extractCP+kk,:)))).^2,2);
                    end
                    vv_temp = mean(metric_d);
                    vv(search_point) = sum(vv_temp(index_k/2+1:end))./sum(vv_temp(1:index_k/2));
                end
                [~, ss] = max(vv);
                timing_offset_est = ss-1;
            end
        otherwise
            error('invalid timing synchronization method.');
    end
    if (subframe_corr == 1)&&(timing_offset_est>0)
        if timing_offset_est+LTE_params.TxSymbols<=length(y_rx)
            y_time_sync = y_rx(timing_offset_est+(1:LTE_params.TxSymbols),:);
        else
            y_time_sync = circshift(y_rx(UE.timing_offset+(1:LTE_params.TxSymbols),:), timing_offset_est - UE.timing_offset);
        end
    else
        y_time_sync = circshift(y_rx, -timing_offset_est);
    end
else
    y_time_sync = y_rx;
    timing_offset_est = 0;
end

if LTE_params.introduce_frequency_offset
    switch UE.freq_sync_method
        case 'perfect'
            if subframe_corr == 1
                y_rx_sync = y_time_sync.*exp(-1i*2*pi*UE.carrier_freq_offset*repmat((timing_offset_est+(0:1:(size(y_time_sync,1)-1)))',1,UE.nRX)/LTE_params.Nfft);
            else
                y_rx_sync = y_time_sync.*exp(-1i*2*pi*UE.carrier_freq_offset*repmat(circshift((0:1:(size(y_time_sync,1)-1))',-timing_offset_est),1,UE.nRX)/LTE_params.Nfft);
            end
            freq_offset_est.error = 0;
            freq_offset_est.int = round(UE.carrier_freq_offset);
            freq_offset_est.frac = UE.carrier_freq_offset - freq_offset_est.int;
            freq_offset_est.res = 0;
        case 'none'
            % for testing the unsynchronized case
            y_rx_sync = y_time_sync;
            freq_offset_est.error = NaN;
            freq_offset_est.int = NaN;
            freq_offset_est.frac = NaN;
            freq_offset_est.res = NaN;
        case 'estimated'
            nTX = size(RefSym,3);
            if (exist('PrimSync','var'))&&(exist('PrimMapping','var'))
                %% estimate the fractional carrier frequency offset
                if (length(LTE_params.Ng) == 2)
                    Index_extractCP_temp{1} = 1:LTE_params.Ng(1);
                    Index_extractCP_temp{2} = reshape(LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
                    Index_extractCP_temp{3} = LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+(1:LTE_params.Ng(1));
                    Index_extractCP_temp{4} = reshape(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+LTE_params.NfftCP{1}+(1:LTE_params.NfftCP{2}*(LTE_params.Ns-1)),LTE_params.NfftCP{2},(LTE_params.Ns-1));
                    Index_extractCP = [Index_extractCP_temp{1} reshape(Index_extractCP_temp{2}(1:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3} reshape(Index_extractCP_temp{4}(1:LTE_params.Ng(2),:),1,[])];
                    %             Index_extractCP = [Index_extractCP_temp{1}(3:end) reshape(Index_extractCP_temp{2}(3:LTE_params.Ng(2),:),1,[]) Index_extractCP_temp{3}(3:end) reshape(Index_extractCP_temp{4}(3:LTE_params.Ng(2),:),1,[])];
                else
                    Index_extractCP_temp = reshape((1:LTE_params.NfftCP*LTE_params.Nsub),LTE_params.NfftCP,LTE_params.Nsub);
                    Index_extractCP = reshape(Index_extractCP_temp(1:LTE_params.Ng,:),1,[]);
                end
                rx_CP1 = y_time_sync(Index_extractCP,:);
                rx_CP2 = y_time_sync(Index_extractCP+LTE_params.Nfft,:);
                expp = sum(mean(rx_CP1 .* conj(rx_CP2)));
                freq_offset_est.frac = -angle(expp)/2/pi;
                if subframe_corr == 1
                    y_rx_sync1 = y_time_sync.*exp(-1i*2*pi*freq_offset_est.frac*repmat((timing_offset_est+(0:1:(size(y_time_sync,1)-1)))',1,UE.nRX)/LTE_params.Nfft);
                else
                    y_rx_sync1 = y_time_sync.*exp(-1i*2*pi*freq_offset_est.frac*repmat(circshift((0:1:(size(y_time_sync,1)-1))',-timing_offset_est),1,UE.nRX)/LTE_params.Nfft);
                end
                %% estimate the integer carrier frequency offset
                PrimSync_rx = zeros(62,UE.nRX);
                SecSync_rx = zeros(62,UE.nRX);
                for nn = 1:UE.nRX
                    y_rx_resolved = cell(4,1);
                    if(length(LTE_params.Ng)==2)
                        y_rx_resolved{1} = y_rx_sync1(1:LTE_params.NfftCP{1},nn);
                        y_rx_resolved{2} = reshape(y_rx_sync1(LTE_params.NfftCP{1}+1:LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
                        y_rx_resolved{3} = y_rx_sync1(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn);
                        y_rx_resolved{4} = reshape(y_rx_sync1(2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:end,nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
                        y_rx_assembled_ifft = [y_rx_resolved{1}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{2}(LTE_params.Index_RxCyclicPrefix{2},:)...
                            y_rx_resolved{3}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{4}(LTE_params.Index_RxCyclicPrefix{2},:)];
                    else
                        y_rx = reshape(y_rx_sync1(:,nn),LTE_params.NfftCP,LTE_params.Nsub);
                        y_rx_assembled_ifft = y_rx(LTE_params.Index_RxCyclicPrefix,:);
                    end
                    y_rx_assembled_shifted = fft(1/sqrt(LTE_params.Nfft)/(sqrt(LTE_params.Nfft/LTE_params.Ntot))*y_rx_assembled_ifft);
                    y_rx_assembled_padded = circshift(y_rx_assembled_shifted,LTE_params.Ntot/2);
                    % remove zero DC carrier
                    y_rx_assembled_temp = y_rx_assembled_padded([1:LTE_params.Ntot/2,LTE_params.Ntot/2+2:LTE_params.Ntot+1],:);
                    PrimSync_rx(:,nn) = y_rx_assembled_temp(PrimMapping);
                    SecSync_rx(:,nn) = y_rx_assembled_temp(SecMapping);
                end
                % ML estimator
                for n_shift = -31:1:31
                    for rr = 1:UE.nRX
                        obs(rr,n_shift+32) = sum(circshift(PrimSync_rx(:,rr).*conj(SecSync_rx(:,rr)),n_shift).*conj(PrimSync).*SecSync);
                    end
                end
                Tp = sum(obs,1).*exp(1i*2*pi*(-31:31)*(LTE_params.NfftCP{2}/LTE_params.Nfft-1));
                [~, argument] = max(Tp);
                freq_offset_est.int = 32 - argument;
                if subframe_corr == 1
                    y_rx_sync = y_rx_sync1.*exp(-1i*2*pi*freq_offset_est.int*repmat((timing_offset_est+(0:1:(size(y_time_sync,1)-1)))',1,UE.nRX)/LTE_params.Nfft);
                else
                    y_rx_sync = y_rx_sync1.*exp(-1i*2*pi*freq_offset_est.int*repmat(circshift((0:1:(size(y_time_sync,1)-1))',-timing_offset_est),1,UE.nRX)/LTE_params.Nfft);
                end
            else
                if subframe_corr == 1
                    y_rx_sync = y_time_sync.*exp(-1i*2*pi*(freq_offset_est.frac+freq_offset_est.int)*repmat((timing_offset_est+(0:1:(size(y_time_sync,1)-1)))',1,UE.nRX)/LTE_params.Nfft);
                else
                    y_rx_sync = y_time_sync.*exp(-1i*2*pi*(freq_offset_est.frac+freq_offset_est.int)*repmat(circshift((0:1:(size(y_time_sync,1)-1))',-timing_offset_est),1,UE.nRX)/LTE_params.Nfft);
                end
            end
            %% estimate the residual carrier frequency offset
            % to frequency domain
            for nn = 1:UE.nRX
                y_rx_resolved = cell(4,1);
                if(length(LTE_params.Ng)==2)
                    y_rx_resolved{1} = y_rx_sync(1:LTE_params.NfftCP{1},nn);
                    y_rx_resolved{2} = reshape(y_rx_sync(LTE_params.NfftCP{1}+1:LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
                    y_rx_resolved{3} = y_rx_sync(LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1),nn);
                    y_rx_resolved{4} = reshape(y_rx_sync(2*LTE_params.NfftCP{1}+LTE_params.NfftCP{2}*(LTE_params.Ns-1)+1:end,nn),LTE_params.NfftCP{2},LTE_params.Ns-1);
                    y_rx_assembled_ifft = [y_rx_resolved{1}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{2}(LTE_params.Index_RxCyclicPrefix{2},:)...
                        y_rx_resolved{3}(LTE_params.Index_RxCyclicPrefix{1},:) y_rx_resolved{4}(LTE_params.Index_RxCyclicPrefix{2},:)];
                else
                    y_rx_block = reshape(y_rx_sync(:,nn),LTE_params.NfftCP,LTE_params.Nsub);
                    y_rx_assembled_ifft = y_rx_block(LTE_params.Index_RxCyclicPrefix,:);
                end
                y_rx_assembled_shifted = fft(1/sqrt(LTE_params.Nfft)/(sqrt(LTE_params.Nfft/LTE_params.Ntot))*y_rx_assembled_ifft);
                y_rx_assembled_padded = circshift(y_rx_assembled_shifted,LTE_params.Ntot/2);
                % remove zero DC carrier
                y_rx_assembled(:,:,nn) = y_rx_assembled_padded([1:LTE_params.Ntot/2,LTE_params.Ntot/2+2:LTE_params.Ntot+1],:);
            end
            % at most reference symbols from 2 TX are used, even if there can be 4 TX
            if nTX~= 4
                nRefUse = nTX;
            else
                nRefUse = 2;
            end
            RefSym = RefSym(:,:,1:nRefUse,:); % extract the reference symbols that are used for RFO estimation
            RefSym_rx = nan([size(RefSym,1), size(RefSym,2), nRefUse, UE.nRX]);   %allocate memory for reference signal for every channel, the number of transmitt antennas is included in RefSym
            for tt = 1:nRefUse
                for rr = 1:UE.nRX
                    RefSym_rx_help = y_rx_assembled(:,:,rr);   %use received symbols from one recieve antenna
                    RefSym_rx_help = RefSym_rx_help(RefMapping(:,:,tt));  %extract the signal on pilots positons
                    RefSym_rx(:,:,tt,rr) = reshape(RefSym_rx_help,size(RefSym(:,:,tt)));  %finally place the reference signal on allocated position
                end
            end
            switch LTE_params.UE_config.rfo_correct_method
                case 'none'
                    freq_offset_est.res = 0;
                case 'subframe'
                    for rr = 1:UE.nRX
                        freq_offset_est.phi(:,:,:,rr) = RefSym_rx(:,[1,2],:,rr).*conj(RefSym_rx(:,[3,4],:,rr)).*conj(RefSym(:,[1,2],:)).*RefSym(:,[3,4],:);
                    end
                    U = sum(sum(sum(sum(freq_offset_est.phi))));
                    freq_offset_est.res = -angle(U)*LTE_params.Nfft/(length(y_rx_sync)*pi);
                otherwise
                    error('unknown RFO estimation method.')
            end
            freq_offset_est.sum = freq_offset_est.frac + freq_offset_est.int + freq_offset_est.res;
            freq_offset_est.error = UE.carrier_freq_offset-freq_offset_est.sum;
            if subframe_corr == 1
                y_rx_sync = y_rx_sync.*exp(-1i*2*pi*freq_offset_est.res*repmat((timing_offset_est+(0:1:(size(y_time_sync,1)-1)))',1,UE.nRX)/LTE_params.Nfft);
            else
                y_rx_sync = y_rx_sync.*exp(-1i*2*pi*freq_offset_est.res*repmat(circshift((0:1:(size(y_time_sync,1)-1))',-timing_offset_est),1,UE.nRX)/LTE_params.Nfft);
            end    
        otherwise
            error('invalid frequency synchronization method.')
    end
else
    y_rx_sync = y_time_sync;
end