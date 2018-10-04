function [LLR_SD,M,H_back,rx_layer_x_equalized] = LTE_detecting(MCS_and_scheduling,BS_nAtPort,rx_user_symbols,UE_nRX,H_est_user,filtering,H_est,LTE_params,receiver,receiver_k,sigma_n2,MSE)
% Detection.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

LLR_SD = 0;
% rx_layer_x = 0;
M = 0;
H_back = 0;
if (~isempty(MCS_and_scheduling.freq_indices))  % check wheter user is assigned RBs
    switch MCS_and_scheduling.tx_mode
        case 1      % SISO
            [LLR_SD,M,rx_layer_x_equalized] = LTE_detect_SIxO(MCS_and_scheduling,filtering,rx_user_symbols,H_est,H_est_user,LTE_params,receiver,receiver_k,sigma_n2,MSE);
        case 2      % Transmit Diversity
            [LLR_SD,M,rx_layer_x_equalized] = LTE_detect_TxD(MCS_and_scheduling,BS_nAtPort,rx_user_symbols,UE_nRX,H_est_user,LTE_params,receiver,receiver_k,sigma_n2);
        case 3      % Open Loop Spatial Multiplexing
            [LLR_SD,M,H_back,rx_layer_x_equalized] = LTE_detect_OLSM(MCS_and_scheduling,rx_user_symbols,H_est_user,LTE_params,BS_nAtPort,receiver,receiver_k,sigma_n2);
        case 4      % Closed Loop Spatial Multiplexing
            [LLR_SD,M,rx_layer_x_equalized] = LTE_detect_CLSM(MCS_and_scheduling,rx_user_symbols,H_est_user,LTE_params,filtering,receiver,receiver_k,sigma_n2);
        case 6      % Interference Alignment
            [LLR_SD,M,H_back,rx_layer_x_equalized] = LTE_detect_IA(MCS_and_scheduling,rx_user_symbols,H_est_user,LTE_params,BS_nAtPort,receiver,receiver_k,sigma_n2);
        otherwise
            error('Mode not supported');
    end
    if LTE_params.UE_config.hard_demapping
        beq_0 = LLR_SD>=0;
        LLR_SD(beq_0)  =  LTE_params.UE_config.LLR_clipping;
        LLR_SD(~beq_0) = -LTE_params.UE_config.LLR_clipping;
    end
end