function [CQI,CQI_bar] = LTE_feedback_SISO(sigma_n2,channel,UE,LTE_params)
% author Stefan Schwarz
% contact stefan.schwarz@nt.tuwien.ac.at
% calculates the CQI feedback for SISO

%% channel prediction
% save channel matrix for channel prediction 
UE.previous_channels = circshift(UE.previous_channels,[0,0,-1,0,0]);
UE.previous_channels(:,:,end,:,:)=channel(:,:,:,:);
H_est_complete = LTE_channel_predictor(UE.previous_channels,LTE_params.uplink_delay,LTE_params.ChanMod_config.filtering,UE.predict);

if LTE_params.feedback.channel_averaging % use channel averaging
    RB_max = LTE_params.Nrb;
    SNR = zeros(RB_max,2); 
else
    RB_max = LTE_params.Ntot;
    SNR = zeros(RB_max,2); 
end

for RB_i = 1:RB_max
    for slot_i = 1:2 
        if LTE_params.feedback.channel_averaging
            freq_band = (RB_i-1)*LTE_params.Nsc+1:min(RB_i*LTE_params.Nsc,size(H_est_complete,1));
        else
            freq_band = RB_i;
        end
%         freq_band = (RB_i-1)*LTE_params.Nsc+1:min(RB_i*LTE_params.Nsc,size(H_est_complete,1));
        H_est = H_est_complete(freq_band,(slot_i-1)*LTE_params.Ns+1:slot_i*LTE_params.Ns,:,:);
        H_t = reshape(mean(mean(H_est,1),2),size(H_est,3),size(H_est,4));   % Channel of current RB and slot
        MSE_mean = mean(UE.MSE(:,freq_band),2);                         % Channel estimator MSE 
        F = pinv(H_t);
        SNR(RB_i,slot_i) = 1./((sigma_n2+MSE_mean(1)).*sum(abs(F).^2,2));
    end    
end
CQI = zeros(LTE_params.Nrb,2);
if LTE_params.UE_config.CQI_fb
    if UE.CQI_fb_gran ~= 1  % single CQI value for whole bandwidth
        CQIs = 0:15;
        SINReff = UE.SINR_averager.average(SNR,CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
        CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
        CQI = CQI_temp*ones(LTE_params.Nrb,2);
        CQI_bar = CQI_temp;
    else
        for RB_ii = 1:LTE_params.Nrb
            if LTE_params.feedback.channel_averaging
                rb_ii = RB_ii;
            else
                rb_ii = (RB_ii-1)*LTE_params.Nsc+1:RB_ii*LTE_params.Nsc;
            end
            for slot_ii = 1:2
                CQIs = (0:15);
                SNR_tmp = SNR(rb_ii,slot_ii);
                SINReff = UE.SINR_averager.average(SNR_tmp(:),CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
    %             SINReff = UE.SINR_averager.average(SNR(RB_ii,slot_ii),CQIs,[LTE_params.CQI_params(CQIs+1).modulation_order]);
                CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
                CQI(RB_ii,slot_ii) = CQI_temp;
            end
        end
        SINReff = UE.SINR_averager.average(SNR,CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
        CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
        CQI_bar = CQI_temp;
    end
else
    CQI = [];
end
