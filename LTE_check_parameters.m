% Check parameter consistency
% Author: Josep Colom, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2009/04/27
% last changes:

if(strcmp(LTE_params.CyclicPrefix,'extended') && LTE_params.SubcarrierSpacing == 7.5e3)
    error('For this combination of subcarrier spacing and cyclic prefix (MBSFN transmissions) reference and synchronization channels not implemented');
end
if LTE_params.max_HARQ_retransmissions > 3 || LTE_params.max_HARQ_retransmissions < 0
    error('Maximum HARQ retransmissions cannot be higher than 3 or negative');
end

if LTE_params.HARQ_processes > 8
    error('The standard does not allow more than 8 HARQ processes. Deleting this error message may result in errors in rate matching process.');
end

% check number of transmit antennas
for bb=1:LTE_params.nBS
    if (BS(bb).nTX == 3 || BS(bb).nTX < 1)
        error('number of antennas not supported');
    end
    if (BS(bb).nAtPort > BS(bb).nTX)
        error('number of antenna ports not consistent with number of antennas');
    end
    for uu = 1:LTE_params.nUE
        if (UE(uu).mode ~= 2 && (LTE_params.scheduler.nLayers(uu) > BS(bb).nTX || LTE_params.scheduler.nLayers(uu) > UE(uu).nRX)) % Exclude TX diversity from this case
            error('number of layers must not be larger than number of antennas');
        end
    end
end

if LTE_params.UE_config.mode == 4
    if LTE_params.BS_config.nTx == 1
        error('Closed loop spatial multiplexing not supported for one transmit antenna');
    end
end    

% Check if the number of subframes is suffucient to estimate the channel
% autocorrelation matrix
% for uu=1:LTE_params.nUE % first assumption: all user equipments have the same capabilities
%     if(UE(uu).realization_num_total >= N_subframes && strcmp(UE(1).autocorrelation_matrix_type,'estimated') && strcmp(UE(1).channel_estimation_method,'MMSE'))
%         error('increase the number of N_subframes or decrease the number of channel realizations used for estimation of channel autocorrelation amtrix')
%     end
% end

if strcmp(LTE_params.scheduler.type,'proportional fair')
    if ~LTE_params.UE_config.CQI_fb
        error('The proportional fair scheduler needs CQI feedback --> activate LTE_params.UE_config.CQI_fb');
    end
    if LTE_params.UE_config.mode == 2 || LTE_params.UE_config.mode == 3
        error('The proportional fair scheduler is not compatible with transmit diversity or open loop spatial multiplexing mode');
    end
    if LTE_params.UE_config.CQI_fb_granularity ~= 1
        disp('A CQI feedback granularity of 1 is recommended for the proportional fair scheduler --> can be set by LTE_params.UE_config.CQI_fb_granularity = 1');
    end
end

if size(SNR_vec,1) < LTE_params.nUE*LTE_params.nBS
    if size(SNR_vec,1) == 1
        SNR_vec = repmat(SNR_vec,LTE_params.nUE*LTE_params.nBS,1);
    else
        error('You have to specify SNR values for all UEs. Every row in the SNR vector corresponds to one UE.');
    end
end

if strcmp(LTE_params.simulation_type,'parallel') && LTE_params.nBS>1
    error('Currently the simulation with multiple eNodeBs is not supported using parellel toolbox');
end
