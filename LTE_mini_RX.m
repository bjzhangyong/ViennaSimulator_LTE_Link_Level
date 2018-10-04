function feedback = LTE_mini_RX(LTE_params, UE_input, SNR, subframe_i, BS_output, UE)
% LTE mini receiver for a specific user.
% [UE_output, ACK, UE] = LTE_RX(chan_output, SNR, AtPort, subframe_i, BS, UE, BS_output)
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   chan_output         ... [1 x 1]struct - channel output:
%           
%                                   [1 x nUE]struct - user specific parameters and variables
% output:   UE_output           ... [1 x nUE]struct UE output for rx_coded_bits bits before decoder and after demapper

%
% date of creation: 2010/11/16


    sigma_n2 = 10^(-SNR/10); %maybe replaced by a function

    
    % Get necessary data and convert in in the correct format
    nLayers    = BS_output.UE_signaling(1).MCS_and_scheduling.nLayers;
    nCodewords = BS_output.UE_signaling(1).MCS_and_scheduling.nCodewords;
    tx_mode    = LTE_params.UE_config.mode;
    

    UE_mapping = ones(LTE_params.Nrb,2); % assuming the user has the full bandwidth
    y_rx_carrier_freq_offset = UE_input.input;
        
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
        rx_user_symbols = zeros(sum(rb_rx_symbols),LTE_params.UE_config.nRX);
        H_est_user = zeros(sum(rb_rx_symbols),LTE_params.UE_config.nRX,LTE_params.BS_config.nTx);
    elseif ~LTE_params.usePDCCH
        rx_symbols  = zeros((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
        PE_SINR_matrix = nan((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
        H_temp_user = zeros((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2,LTE_params.UE_config.nRX,LTE_params.BS_config.nTx);
    else
        rb_rx_symbols = LTE_params.Nsc*LTE_params.Ns - total_no_refsym -  CHusedElements(UE_mapping);
        rx_user_symbols = zeros(sum(rb_rx_symbols),LTE_params.UE_config.nRX);
        H_est_user = zeros(sum(rb_rx_symbols),LTE_params.UE_config.nRX,LTE_params.BS_config.nTx);
    end
    
%     %% Carrier Frequency Offset Compensation
%     if LTE_params.introduce_frequency_offset
%         if (subframe_corr == 1 || subframe_corr == 6)
%             [y_rx_sync, UE_output(uu).freq_offset_est] = LTE_rx_freq_sync(LTE_params, y_rx_carrier_freq_offset, UE(uu), UE_output(uu).freq_offset_est, RefSym, RefMapping, PrimSync, PrimMapping, SecSync, SecMapping);
%         else
%             [y_rx_sync, UE_output(uu).freq_offset_est] = LTE_rx_freq_sync(LTE_params, y_rx_carrier_freq_offset, UE(uu), UE_output(uu).freq_offset_est, RefSym, RefMapping);
%         end
%     else
%         y_rx_sync = UE_input.input;
%         UE_output(uu).freq_offset_est.error = NaN;
%         UE_output(uu).freq_offset_est.frac = NaN;
%         UE_output(uu).freq_offset_est.int = NaN;
%         UE_output(uu).freq_offset_est.res = NaN;
%     end
    y_rx_sync = UE_input.input;
    
    %% Remove CP, FFT, remove zeros
    for nn = 1:LTE_params.UE_config.nRX
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
    rx_ref_symbols = nan([size(RefSym),LTE_params.UE_config.nRX]);   %allocate memery for reference signal for every channel, the number of transmitt antennas is included in RefSym
    for tt = 1:LTE_params.BS_config.nTx
        for rr = 1:LTE_params.UE_config.nRX
            rx_ref_symbols_help = y_rx_assembled(:,:,rr);   %use recieved symbols from one recieve antenna
            rx_ref_symbols_help = rx_ref_symbols_help(RefMapping(:,:,tt));  %extract the signal on pilots positons
            if(tt>2)    %in case of 3rd and 4th transmitt antenna, there are less pilots, so we have to insert some nan symbols,to be consistent
                rx_ref_symbols_help = [rx_ref_symbols_help;nan(size(rx_ref_symbols_help))];
            end
            rx_ref_symbols(:,:,tt,rr) = reshape(rx_ref_symbols_help,size(RefSym(:,:,tt)));  %finally place the reference signal on allocated position
        end
    end
    
    %% Channel and noise estimation
    H_est = LTE_channel_estimator_mini(LTE_params, rx_ref_symbols, RefSym, RefMapping, sigma_n2, y_rx_assembled);
    

    
%% calculation of feedback
[RI_tmp,PMI_tmp,CQI_tmp, CQI_bar]=LTE_feedback(LTE_params.BS_config.nTx,sigma_n2,LTE_params,H_est,UE(1),1,LTE_params.UE_config.mode);
feedback = 0;
    
% 
%     % Feedback performed only when data is received
%     feedback_rv_idx = [UE_signaling.turbo_rate_matcher.rv_idx];
%     UE_output(uu).rv_idx = feedback_rv_idx;
%     
%     % Precoding feedback calculation
%     if LTE_params.uplink_delay ~= 0
%         %             H_est = LTE_channel_estimator(LTE_params,ChanMod, UE(uu), rx_ref_symbols, RefSym, RefMapping, ChanMod_output.genie.H_fft,sigma_n2,BS_output.cell_genie.y_tx_assembled,y_rx_assembled);
%         %             [RI_tmp,PMI_tmp]=LTE_feedback_precoding(BS.nAtPort,sigma_n2,LTE_params,UE_signaling.MCS_and_scheduling(1).CQI_params(1).modulation_order(1),ChanMod_output.genie.H_fft,UE(uu),uu);
%         %             if LTE_params.UE_config.mode == 3
%         %                 [RI_tmp,PMI_tmp,CQI_tmp]=LTE_feedback(BS.nAtPort,sigma_n2,LTE_params,H_back,UE(uu),uu,LTE_params.UE_config.mode);
%         %             else
%         [RI_tmp,PMI_tmp,CQI_tmp]=LTE_feedback(BS.nAtPort,sigma_n2,LTE_params,H_est,UE(uu),uu,LTE_params.UE_config.mode);
%         %             end
%         UE_output(uu).RI = RI_tmp;
%         UE_output(uu).PMI = PMI_tmp-1;
%         UE_output(uu).CQI = CQI_tmp;
%     end
% 
%     
% 
% 
% % Finish feedback construction
% % For now use perfect channel information. Can be changed to an SINR computed from channel estimation
% % SINR = BS_output.cell_genie.SINR(uu,:);
% % SINRs_to_CQI = SINR(1:12:end);
% % UE_output(uu).CQI_feedback = floor(LTE_common_CQI_mapping(LTE_params.CQI_mapping,SINRs_to_CQI)); % only one CQI value over bandwidth, because there is only noise added
% % UE_output(uu).CQI = UE_output(uu).CQI_feedback;
% % UE_output(uu).PMI = 1;
% % HARQ process id of this current ACK
% UE_output(uu).HARQ_process = BS_output.UE_signaling(uu).MCS_and_scheduling.HARQ_process_id;
% UE_output(uu).UE_scheduled = (BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs>0);
% 
% %% Calculate noise on a subcarrier basis with the genie information
% % noise_SC = BS_output.cell_genie.y_tx_assembled - y_rx_assembled;
% % SINR_SC_linear = mean(abs(noise_SC).^2,2);
% % BS_output.cell_genie.genie_SNR = SINR_SC_linear;
