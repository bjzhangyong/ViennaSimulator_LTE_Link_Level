% generate UE, BS ...
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2010/01/19
% last changes: 

%% Generate UEs
UE = network_elements.UE(LTE_params.UE_config,LTE_params.Ntot,LTE_params.Nsub,LTE_params.BS_config.nTx,LTE_params.ChanMod_config.filtering);
for uu=1:LTE_params.nUE*LTE_params.nBS % first assumption: all user equipments have the same capabilities
    UE(uu) = network_elements.UE(LTE_params.UE_config,LTE_params.Ntot,LTE_params.Nsub,LTE_params.BS_config.nTx,LTE_params.ChanMod_config.filtering);
end

%% Append the traffic model to the UEs
if LTE_params.trafficmodel.usetraffic_model && ~strcmp(LTE_params.scheduler.type,'constrained scheduler')
    warning('Traffic models are just utilized in connection with the constrained scheduler - deactivating traffic models!');
    LTE_params.trafficmodel.usetraffic_model = false;
end
for uu = 1:LTE_params.nUE*LTE_params.nBS
    if exist('model','var')
        if ~isempty(model)
            UE(uu).traffic_model = LTE_trafficmodel(LTE_params.trafficmodel,UE(uu),max(1,LTE_params.uplink_delay),model(uu));
        end
    else
        UE(uu).traffic_model = LTE_trafficmodel(LTE_params.trafficmodel,UE(uu),max(1,LTE_params.uplink_delay));
    end
end

%% Generate eNodeBs
LTE_params.BS_config.AtPort = 0:LTE_params.BS_config.nTx-1;
LTE_params.BS_config.nAtPort = length(LTE_params.BS_config.AtPort);
LTE_params.BS_config.maxStreams = 2;
LTE_params.BS_config.HARQ_processes = LTE_params.HARQ_processes;

for bb=1:LTE_params.nBS % first assumption: all base stations have the same capabilities
    BS(bb) = network_elements.eNodeB(...
        bb,...
        LTE_params.BS_config.nTx,...
        LTE_params.nUE,...
        LTE_params.BS_config.AtPort,...
        LTE_params.BS_config.nAtPort,...
        LTE_params.BS_config.maxStreams,...
        LTE_params.BS_config.HARQ_processes,...
        LTE_params.max_HARQ_retransmissions);
end

%% Pregenerate Reference and Synchronization signals for a whole frame

for b_ = 1:length(BS)
    for subframe_num = 1:10
        % Reference signals
        [LTE_params.Reference_Signal(b_,subframe_num).RefSym,...
            LTE_params.Reference_Signal(b_,subframe_num).RefMapping,...
            LTE_params.Reference_Signal(b_,subframe_num).total_no_refsym,...
            LTE_params.Reference_Signal(b_,subframe_num).NoData_indices] = LTE_common_gen_Reference_Signal(LTE_params.Nrb, LTE_params.Nsub, BS(b_).nTX, BS(b_).NIDcell, subframe_num);
        % Synchronization signals
        [LTE_params.Sync_Signal(b_,subframe_num).PrimSync,...
            LTE_params.Sync_Signal(b_,subframe_num).PrimMapping,...
            LTE_params.Sync_Signal(b_,subframe_num).SecSync,...
            LTE_params.Sync_Signal(b_,subframe_num).SecMapping,...
            LTE_params.Sync_Signal(b_,subframe_num).SyncUsedElements,...
            LTE_params.Sync_Signal(b_,subframe_num).ResMapping,...
            LTE_params.Sync_Signal(b_,subframe_num).ResUsedElements,...
            LTE_params.Sync_Signal(b_,subframe_num).ResNElements] = LTE_common_gen_Synchronization_Signal(BS(b_).NIDcell, LTE_params.Nrb, LTE_params.Nsc, LTE_params.Nsub, subframe_num);
    
        % Pilot symbol boosting
        if ~LTE_params.Pilots_power_offset==0 
            total_number_of_pilots = sum(sum(sum(LTE_params.Reference_Signal(1,subframe_num).RefMapping)));
            total_number_of_data = LTE_params.Ntot*LTE_params.Nsub - total_number_of_pilots;
            LTE_params.pilot_factor = (LTE_params.Ntot*LTE_params.Nsub)/(10^(-LTE_params.Pilots_power_offset/20)*total_number_of_data + total_number_of_pilots);
            LTE_params.data_factor = 10^(-LTE_params.Pilots_power_offset/20) * LTE_params.pilot_factor;

            LTE_params.Reference_Signal(1,subframe_num).RefSym = LTE_params.pilot_factor*LTE_params.Reference_Signal(1,subframe_num).RefSym;                    
        end
    end
end

%% Pregenerate other channels (PBCH,PDCCH)
for b_ = 1:length(BS)
    for subframe_num = 1:10    
         [LTE_params.PBCH(b_,subframe_num).Mapping,LTE_params.PBCH(b_,subframe_num).UsedElements,LTE_params.PBCH(b_,subframe_num).N_Elements,NoData] = LTE_common_gen_PBCH(LTE_params.Nrb,LTE_params.Nsc,LTE_params.Nsub,BS(b_).NIDcell,subframe_num,LTE_params.Ns,BS(b_).nTX);
         [LTE_params.PDCCH(b_,subframe_num).Mapping,LTE_params.PDCCH(b_,subframe_num).UsedElements,LTE_params.PDCCH(b_,subframe_num).N_Elements] = LTE_common_gen_PDCCH(LTE_params.Nrb,LTE_params.Nsc,LTE_params.Nsub,LTE_params.Ns,LTE_params.Bandwidth,NoData);   
    end
end
LTE_params.scheduler.overhead_ref = LTE_params.Reference_Signal(1).total_no_refsym;
LTE_params.scheduler.overhead_sync = mean((LTE_params.Sync_Signal(1).ResUsedElements+LTE_params.Sync_Signal(1).SyncUsedElements(:)));
if LTE_params.usePBCH 
    LTE_params.scheduler.overhead_sync = LTE_params.scheduler.overhead_sync + mean([LTE_params.PBCH.N_Elements]/LTE_params.Nrb);
end
if LTE_params.usePDCCH
    LTE_params.scheduler.overhead_sync = LTE_params.scheduler.overhead_sync+mean([LTE_params.PDCCH.N_Elements]/LTE_params.Nrb);
end

%% update number of channel realizations
for uu=1:LTE_params.nUE*LTE_params.nBS
    UE(uu).realization_num_total = UE(uu).realization_num_total * UE(uu).nRX * ChanMod.nTX; %channel of every subframe consitcs of nRx*nTx different channel
end

%% Load or allocate channel autocorrelation matrix
for uu=1:LTE_params.nUE*LTE_params.nBS
     if strcmp(UE(uu).channel_estimation_method,'MMSE') || strcmp(UE(uu).channel_estimation_method,'ALMMSE2') || strcmp(UE(uu).channel_estimation_method,'ALMMSE') || strcmp(UE(uu).channel_estimation_method,'LS')
%         if strcmp(UE(uu).autocorrelation_matrix_type,'ideal')
            %load autocorrelation_matrix_PedB;  
            switch ChanMod.interpolation_method
                case 'shift_to_nearest_neighbor'
                    switch LTE_params.ChanMod_config.type
                        case {'AWGN','flat Rayleigh'}
                              H = 1; 
                        otherwise
                            NTap = size(ChanMod.PDP_dB,2);% Number of Taps
                            tap_delays = round(ChanMod.PDP_dB(2,:)*LTE_params.Fs);
                            H = zeros(tap_delays(end)+1,1);
            %                     for tap_i = 1 : NTap
            %                         if(tap_i == 1)
            %                             prev_delay = -1;
            %                         else
            %                             prev_delay = tap_delays(tap_i-1);
            %                         end
            %                         tap_delays(tap_i)
            %                         if(tap_delays(tap_i) ~= prev_delay) % when some taps merge by a specific sampling freq
            %                             H(tap_delays(tap_i)+1) = sqrt(10.^(ChanMod.PDP_dB(1,tap_i)./10))*1;
            %                         end
            %                     end
                            for tap_i = 1:tap_delays(end)+1
                                H(tap_i) = sqrt(sum(10.^(ChanMod.PDP_dB(1,tap_delays == tap_i-1)./10)));
                            end
                            H = H./ChanMod.normH;    
                    end
                case 'sinc_interpolation'
                    H = sum(ChanMod.interpolator.precomputed_sincs)./ChanMod.interpolator.norm_h;
            end
            DFT_matrix = dftmtx(LTE_params.Nfft);
            time_correlation_matrix = zeros(LTE_params.Nfft);
            time_correlation_matrix(1:length(H),1:length(H)) = diag(H.^2);
            % calculate perfect autocorrelation function R_ff = DFT_mat * R_tt * DFT_mat'
            autocorrelation_matrix_help = DFT_matrix * time_correlation_matrix * DFT_matrix';          
            % pick the correct values, remove DC ...
%             [LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]
%             autocorrelation_matrix = autocorrelation_matrix_help([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 1:LTE_params.Ntot/2],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 1:LTE_params.Ntot/2]);
            autocorrelation_matrix = autocorrelation_matrix_help([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
            UE(uu).channel_autocorrelation_matrix = autocorrelation_matrix;
        if strcmp(UE(uu).autocorrelation_matrix_type,'estimated')
            UE(uu).channel_autocorrelation_matrix= zeros(LTE_params.Ntot);
%         else
%             error('wrong type of autocorrelation matrix');
        end
    end
end

%% Channel parameters dependent - now only the same channel parameters for each user and BS are allowed
% load Correlation Matrices
if strcmp(ChanMod.type,'PedA') || strcmp(ChanMod.type,'PedB') || strcmp(ChanMod.type,'VehA') || strcmp(ChanMod.type,'VehB') ||strcmp(ChanMod.type,'TU') || strcmp(ChanMod.type,'RA') || strcmp(ChanMod.type,'HT') || strcmp(ChanMod.type,'Rayleigh2')
    ChanMod.corrRX = ones(size(ChanMod.PDP_dB,2),UE(1).nRX,UE(1).nRX);
    ChanMod.corrTX = ones(size(ChanMod.PDP_dB,2),BS(1).nTX,BS(1).nTX);
    for kk = 1:size(ChanMod.PDP_dB,2)
        ChanMod.corrRX(kk,:,:) = eye(UE(1).nRX);
        ChanMod.corrTX(kk,:,:) = eye(BS(1).nTX);
    end
elseif strcmp(ChanMod.type,'PedBcorr') || strcmp(ChanMod.type,'EVehA') || strcmp(ChanMod.type,'ETU')
    ChanMod.corrRX = ones(size(ChanMod.PDP_dB,2),UE(1).nRX,UE(1).nRX);
    ChanMod.corrTX = ones(size(ChanMod.PDP_dB,2),BS(1).nTX,BS(1).nTX);
    for kk = 1:size(ChanMod.PDP_dB,2)
        ChanMod.corrRX(kk,:,:) = eye(UE(1).nRX) + ChanMod.corr_coefRX*ones(UE(1).nRX) - ChanMod.corr_coefRX*eye(UE(1).nRX);
        ChanMod.corrTX(kk,:,:) = eye(BS(1).nTX) + ChanMod.corr_coefTX*ones(BS(1).nTX) - ChanMod.corr_coefTX*eye(BS(1).nTX);
    end
end

%% generate psi, theta, phi for Rosa-Zheng Channel Model
switch ChanMod.time_correlation 
    case 'correlated'
    switch ChanMod.type
        case {'PedA', 'PedB', 'PedBcorr','VehA','VehB','TU','RA','HT','EPedA','EVehA','ETU'}
            for uu=1:LTE_params.nUE*LTE_params.nBS

                number_of_taps = ChanMod.interpolator.num_faders;

                UE(uu).channel_coef_rosa_zheng.theta = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps)*2 -1) * pi;
                UE(uu).channel_coef_rosa_zheng.phi   = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi; % The Original paper states \phi as not changing for each sinusoid, but that is not the case (see T. Zemen and C.F. Mecklenbr¨auker, “Time-Variant Channel Estimation Using Discrete Prolate Spheroidal Sequences,” IEEE Transactions on Signal Processing, vol. 53, no. 9, pp. 3597–3607, Sept. 2005.)
                UE(uu).channel_coef_rosa_zheng.psi   = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;

            end
        otherwise
    end
    case 'independent'
end

%% Scheduler initial initialization

if strcmp(LTE_params.scheduler.assignment,'static')
    LTE_params.scheduler.PMI = ones(LTE_params.Nrb,2)*LTE_params.scheduler.PMI;
else
    LTE_params.scheduler.PMI = ones(LTE_params.Nrb,2)*LTE_params.scheduler.PMI;
end

% if(strcmp(LTE_params.scheduler.cqi,'set') && strcmp(LTE_params.scheduler.assignment,'static') && ~strcmp(LTE_params.scheduler.type,'best cqi')) || strcmp(LTE_params.scheduler.assignment,'semi static') || strcmp(LTE_params.scheduler.type,'proportional fair')
if(strcmp(LTE_params.scheduler.cqi,'set') && strcmp(LTE_params.scheduler.assignment,'static'))  || strcmp(LTE_params.scheduler.assignment,'semi static') || strcmp(LTE_params.scheduler.type,'proportional fair') || strcmp(LTE_params.scheduler.type,'best cqi') || strcmp(LTE_params.scheduler.type,'max throughput')
    LTE_params.scheduler.cqi = cqi_i;
elseif strcmp(LTE_params.scheduler.type,'fixed')
    LTE_params.scheduler.cqi = cqi_i;
    LTE_params.scheduler.assignment = 'static';
else
    LTE_params.scheduler.assignment = 'dynamic';
end

if LTE_params.uplink_delay==0
    LTE_params.scheduler.zero_delay = true;
else
    LTE_params.scheduler.zero_delay = false;
end

LTE_params.scheduler.UE_specific = BS(1).UE_specific;

% Copy CQI mapping configuration to the scheduler to allow it to call the CQI mapping funciton with correct parameters
LTE_params.scheduler.CQI_mapping_params = LTE_params.CQI_mapping;

LTE_params.scheduler.CDD = LTE_params.UE_config.CDD;
for uu=1:LTE_params.nUE*LTE_params.nBS
    for i =1:length(LTE_params.nUE)
        LTE_params.scheduler.tx_mode(uu) = LTE_params.UE_config.mode;      % assigned trasmission mode for each user
        switch UE(uu).mode
            case 1 % Single Antenna
                LTE_params.scheduler.nLayers(uu)  = 1;
                LTE_params.scheduler.nCodewords(uu) = 1;
            case 2 % Transmit Diversity
                LTE_params.scheduler.nLayers(uu)  = BS.nAtPort;
                LTE_params.scheduler.nCodewords(uu) = 1;
            case 3 % Open Loop Spatial Multiplexing
                if exist('RANlayers','var')
                    LTE_params.scheduler.nLayers(uu)  = RANlayers;
                else
                    LTE_params.scheduler.nLayers(uu)  = min([LTE_params.BS_config.nTx LTE_params.UE_config.nRX]);
                end
                LTE_params.scheduler.nCodewords(uu) = min(2,LTE_params.scheduler.nLayers(uu));
            case 4 % Closed Loop SM
                LTE_params.scheduler.nLayers(uu)  = min([LTE_params.BS_config.nTx LTE_params.UE_config.nRX]); 
%                 LTE_params.scheduler.nLayers(uu) =  3;
                LTE_params.scheduler.nCodewords(uu) =  min(2,LTE_params.scheduler.nLayers(uu));
            case 6 % Interference Alignment
                LTE_params.scheduler.nLayers(uu) = LTE_params.IA_streams(uu);
                LTE_params.scheduler.nCodewords(uu) = 1;
            otherwise
                error('TX mode not supported. TX mode %d found',UE(uu).mode);
        end
    end
end

%% Further scheduler initialization: Create scheduler. Assume the eNodeB where the scheduler is to be put is the first one
LTE_params.scheduler.max_HARQ_retx = LTE_params.max_HARQ_retransmissions;

for bb = 1:LTE_params.nBS
    switch LTE_params.scheduler.type
        case 'round robin'
            BS(bb).scheduler = network_elements.roundRobinScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'fixed'
            averager = network_elements.eesmAverager(LTE_params.UE_config.SINR_averaging.EESMbetas,0:15);
            if strcmp(LTE_params.scheduler.assignment,'semi static')
                mapping_data = LTE_params.CQI_mapping;
                BS(bb).scheduler = network_elements.AdaptiveFeedbackScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data);
            else
                mapping_data = LTE_params.CQI_mapping.table;
                BS(bb).scheduler = network_elements.fixedScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data);
            end
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'best cqi'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = network_elements.bestCqiSchedulerMIMO(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'proportional fair'
            % %         IPF
            %         averager = network_elements.eesmAverager(LTE_params.UE_config.SINR_averaging.EESMbetas,0:15);
            %         mapping_data = LTE_params.CQI_mapping.table;
            %         BS(1).scheduler = network_elements.InstProportionalFairScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE,LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data);
            
            % %         PF
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.PropFair_Sun(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order],LTE_params.Tsubframe);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'max throughput'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.MaxThroughputScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'max min'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.MaxMinScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'resource fair'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.ResourceFairScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'Kwan'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.Kwan_scheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'opt throughput'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.OptimumThroughputScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'var fair'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.VariableFairnessScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'alpha fair'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.PropFair_Sun_Var(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order],LTE_params.Tsubframe);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'weighted sum rate'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.Weighted_sumrate(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order],LTE_params.scheduler.weights);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'variable fair'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.VarFairScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'best cqi HARQ single-user'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = network_elements.bestCqiSchedulerMIMOHARQ(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'utility max'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.UtilMaxScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        case 'constrained scheduler'
            averager = network_elements.miesmAverager(LTE_params.UE_config.SINR_averaging.MIESMbetas,0:15);
            mapping_data = LTE_params.CQI_mapping;
            BS(bb).scheduler = schedulers.ConstrainedScheduler(LTE_params.Nrb,LTE_params.Nsc*LTE_params.Ns,UE(LTE_params.connection_table(bb,:)),LTE_params.scheduler,LTE_params.CQI_params,averager,mapping_data,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(1:15).modulation_order]);
            BS(bb).scheduler.attached_eNodeB = BS(bb);
        otherwise
            error('scheduler not supported');
    end
end

MI_data_loader = utils.MIdata_loader; 
LTE_params.MI_data = MI_data_loader.MIdata_load;
clear MI_data_loader
for i = 1:3
    LTE_params.MI_data(i).BICM = interp1(LTE_params.MI_data(i).SNR,LTE_params.MI_data(i).BICM,linspace(LTE_params.MI_data(i).SNR(1),LTE_params.MI_data(i).SNR(end),1000));
    LTE_params.MI_data(i).SNR = linspace(LTE_params.MI_data(i).SNR(1),LTE_params.MI_data(i).SNR(end),1000);
end
