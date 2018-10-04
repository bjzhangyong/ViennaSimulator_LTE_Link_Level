function [rank_i,PMI,CQI,CQI_bar] =  LTE_feedback(nAtPort,sigma_n2,LTE_params,channel,UE,uu,modus,cqi_i)
% author Stefan Schwarz
% contact stefan.schwarz@nt.tuwien.ac.at
% calculates the PMI, RI and CQI feedback

switch modus
    case 1
        if UE.CQI_fb 
            [CQI,CQI_bar] = LTE_feedback_SISO(sigma_n2,channel,UE,LTE_params);
        else
            CQI = [];
            CQI_bar = [];
        end
        rank_i = [];
        PMI = [];
    case 2
%         rank_i = [];
%         PMI = [];
%         CQI = cqi_i*ones(LTE_params.Nrb,2);
%         CQI_bar = 0;
        [CQI] = LTE_feedback_TxD(nAtPort,sigma_n2,LTE_params,channel,UE,uu);
        CQI_bar = 0;
        PMI = [];
        rank_i = [];
    case {3,6}
%         if ~LTE_params.uplink_delay
%             channel = reshape(channel,size(channel,1)*size(channel,2),size(channel,3),size(channel,4));
%         end
%         [CQI] = LTE_feedback_OLSM(sigma_n2,channel,UE,LTE_params);
%         rank_i = [];
%         PMI = [];
%         CQI_bar = 0;
        [rank_i,CQI] = LTE_feedback_OLSM(nAtPort,sigma_n2,LTE_params,channel,UE,uu);
        CQI_bar = 0;
        PMI = 0;
    case 4
        [rank_i,PMI,CQI] = LTE_feedback_CLSM(nAtPort,sigma_n2,LTE_params,channel,UE,uu);
%         [rank_i,PMI,CQI] = LTE_feedback_CLSM_test(nAtPort,sigma_n2,LTE_params,channel,UE,uu);
        CQI_bar = 0;
end

