function [CQI] =  LTE_feedback_TxD(nAtPort,sigma_n2,LTE_params,channel,UE,uu)
% author Stefan Schwarz
% contact stefan.schwarz@nt.tuwien.ac.at
% calculates the PMI, RI and CQI feedback

%% channel prediction
% save channel matrix for channel prediction 
UE.previous_channels = circshift(UE.previous_channels,[0,0,-1,0,0]);
UE.previous_channels(:,:,end,:,:)=channel(:,:,:,:);
H_est_complete = LTE_channel_predictor(UE.previous_channels,LTE_params.uplink_delay,LTE_params.ChanMod_config.filtering,UE.predict);

%% some initializations
% rank_i = min(size(H_est_complete,3),size(H_est_complete,4));    % NOTE: remove when rank_i feedback is supported!
rank_loop = LTE_params.scheduler.nLayers(uu);
H = reshape(mean(mean(H_est_complete,1),2),size(H_est_complete,3),size(H_est_complete,4)); 
CQI = zeros(LTE_params.Nrb,2,rank_loop);
CQI_temp = zeros(LTE_params.Nrb,2,rank_loop);

% if nAtPort == 2 
%     i_max = [3,2];
% else
%     i_max = [15,15,15,15];
% end

if LTE_params.feedback.channel_averaging % use channel averaging
    RB_max = LTE_params.Nrb;
    I_g = zeros(RB_max,2);   % total sum rate over all resource blocks 
    SNR_g = zeros(RB_max,2,4); 
else
    RB_max = LTE_params.Ntot;
    I_g = zeros(RB_max,2);   % total sum rate over all resource blocks 
    SNR_g = zeros(RB_max,2,4); 
end
% figure(1)
% clf
for RB_i = 1:RB_max
    for slot_i = 1:2  
        if LTE_params.feedback.channel_averaging
            freq_band = (RB_i-1)*LTE_params.Nsc+1:min(RB_i*LTE_params.Nsc,size(H_est_complete,1));
        else
            freq_band = RB_i;
        end
        H_est = H_est_complete(freq_band,(slot_i-1)*LTE_params.Ns+1:slot_i*LTE_params.Ns,:,:);
        
                H_t = reshape(mean(mean(H_est,1),2),size(H_est,3),size(H_est,4));   % Channel of current RB and slot 
        if ~LTE_params.feedback.ignore_channel_estimation % channel estimation error will be ignored if set so
            MSE_mean = mean(UE.MSE(:,freq_band),2);                         % Channel estimator MSE 
        else
            MSE_mean = zeros(size(UE.MSE,1),1);
        end
%         I = zeros(1,1);   
        SNR = zeros(max(rank_loop),1);
        for rr = rank_loop    % Iterate over all possible layer numbers, to find the one that delivers the maximum sum rate   
%             if rr == 2 && nAtPort == 2 
%                 i_start = 1;
%                 i_stop = i_max(rr);
%             else
%                 i_start = 0;
%                 i_stop = i_max(rr);
%             end
%             for i = i_start:i_stop
%                 [Z W U D] = LTE_common_get_precoding_matrix(2,nAtPort,0,rr,LTE_params);
%                 P = H_t*WDU;

                switch size(H,1)
                    case 1
                        switch size(H,2)
                            case 2
                                P = [H_t(1,1),-H_t(1,2);conj(H_t(1,2)),conj(H_t(1,1))];
                            case 4
                                P = [H_t(1,1),-H_t(1,3),0,0;conj(H_t(1,3)),conj(H_t(1,1)),0,0;...
                                    0,0,H_t(1,2),-H_t(1,4);0,0,conj(H_t(1,4)),conj(H_t(1,2))];
                            otherwise
                                error('not implemented')
                        end
                    case 2
                        switch size(H,2)
                            case 2
                                P = [H_t(1,1),-H_t(1,2);H_t(2,1),-H_t(2,2);conj(H_t(1,2)),conj(H_t(1,1));conj(H_t(2,2)),conj(H_t(2,1))]; 
                            case 4
                                P = [H_t(1,1),-H_t(1,3),0,0;H_t(2,1),-H_t(2,3),0,0;conj(H_t(1,3)),conj(H_t(1,1)),0,0;conj(H_t(2,3)),conj(H_t(2,1)),0,0;...  % generate equivalent channel matrix
                                    0,0,H_t(1,2),-H_t(1,4);0,0,H_t(2,2),-H_t(2,4);0,0,conj(H_t(1,4)),conj(H_t(1,2));0,0,conj(H_t(2,4)),conj(H_t(2,2))];
                            otherwise
                                error('not implemented')
                        end
                    case 4
                        switch size(H,2)
                            case 4
%                                P = [H_t(1,1),-H_t(1,3),0,0;H_t(2,1),-H_t(2,3),0,0;conj(H_t(1,3)),conj(H_t(1,1)),0,0;conj(H_t(2,3)),conj(H_t(2,1)),0,0;...  % generate equivalent channel matrix
%                                     0,0,H_t(1,2),-H_t(1,4);0,0,H_t(2,2),-H_t(2,4);0,0,conj(H_t(1,4)),conj(H_t(1,2));0,0,conj(H_t(2,4)),conj(H_t(2,2));...
%                                     H_t(3,1),-H_t(3,3),0,0;H_t(4,1),-H_t(4,3),0,0;conj(H_t(3,3)),conj(H_t(3,1)),0,0;conj(H_t(4,3)),conj(H_t(4,1)),0,0;...  % generate equivalent channel matrix
%                                     0,0,H_t(3,2),-H_t(3,4);0,0,H_t(4,2),-H_t(4,4);0,0,conj(H_t(3,4)),conj(H_t(3,2));0,0,conj(H_t(4,4)),conj(H_t(4,2))];
                                P = [H_t(1,1),-H_t(1,3),0,0;H_t(2,1),-H_t(2,3),0,0;H_t(3,1),-H_t(3,3),0,0;H_t(4,1),-H_t(4,3),0,0;conj(H_t(1,3)),conj(H_t(1,1)),0,0;conj(H_t(2,3)),conj(H_t(2,1)),0,0;conj(H_t(3,3)),conj(H_t(3,1)),0,0;conj(H_t(4,3)),conj(H_t(4,1)),0,0;...  % generate equivalent channel matrix
                                    0,0,H_t(1,2),-H_t(1,4);0,0,H_t(2,2),-H_t(2,4);0,0,H_t(3,2),-H_t(3,4);0,0,H_t(4,2),-H_t(4,4);0,0,conj(H_t(1,4)),conj(H_t(1,2));0,0,conj(H_t(2,4)),conj(H_t(2,2));0,0,conj(H_t(3,4)),conj(H_t(3,2));0,0,conj(H_t(4,4)),conj(H_t(4,2))];
                            otherwise
                                error('not implemented')
                        end 
                    otherwise
                        error('not implemented')
                end
                P = P*1/sqrt(2);
                switch LTE_params.UE_config.receiver 
%                   switch 'ZF'
                    case 'ZF'
                        F = pinv(P);   % ZF receiver                  
                    case 'MMSE'
                        temp = P'*P;
                        F = (temp+sigma_n2*eye(size(temp)))^-1*P';  % MMSE receiver
%                     case 'SSD'
%                         Jr = real(H_t*W);
%                         Ji = imag(H_t*W);
%                         J = [Jr(1,1),-Ji(1,1),Jr(1,2),-Ji(1,2);Ji(1,1),Jr(1,1),Ji(1,2),Jr(1,2);Jr(2,1),-Ji(2,1),Jr(2,2),-Ji(2,2);Ji(2,1),Jr(2,1),Ji(2,2),Jr(2,2)]; % Lattice generator
%                         % calculate error covariance matrix, assuming AWGN
%                         E(:,:,i) = (J'*J)^-1;
%                         F = eye(size(P,1));
                    otherwise
                        temp = P'*P;
                        F = (temp+sigma_n2*eye(size(temp)))^-1*P';  % MMSE receiver
                end
                K = F*P;
                SNR(:) = abs(diag(K)).^2./(sum(abs(K-diag(diag(K))).^2,2)+(sigma_n2+MSE_mean(1:rr)).*sum(abs(F).^2,2));
%                 I(i+1,rr) = sum(log2(1+SNR(i+1,rr,:)));  % rate of one resource block for different precoders and ranks
%                 I(rr) = LTE_feedback_getBICM(LTE_params,squeeze(SNR(rr,:)));
%             end
        end
%         I_g(RB_i,slot_i) = I_g(RB_i,slot_i) + I;
        SNR = 10*log10(SNR);
        SNR(isnan(SNR))= -inf;
        SNR_g(RB_i,slot_i,1:max(rank_loop)) = SNR;
%         if UE.PMI_fb_gran == 1 
%             [~,C1] = max(I,[],1);
%             PMI_temp(:,RB_i,slot_i) = C1;
%         end
    end
end

% I_ges = sum(sum(I_g,2),3);
% [~,I1] = max(max(I_ges,[],1));
rank_i = rank_loop;    % choose rank indicator to maximize mutual information
if LTE_params.UE_config.CQI_fb % wheter CQI feedback is used
%     for i2 = 1:min(2,rank_i)   % iterate over the number of streams
        if UE.CQI_fb_gran ~= 1  % single CQI value for whole bandwidth
%             if rank_i == 4
%                 AWGN_SNR = squeeze(SNR_g(:,:,(i2-1)*2+1:i2*2));
%             elseif rank_i == 3 && i2 == 2
%                 AWGN_SNR = squeeze(SNR_g(:,:,i2:rank_i));
%             else
                AWGN_SNR = squeeze(SNR_g(:,:,1:rank_i));
%             end
            CQIs = (0:15);
%                 SINReff = UE.SINR_averager.average(10.^(AWGN_SNR(:)/10),CQIs,[LTE_params.CQI_params(CQIs+1).modulation_order]);
            SINReff = UE.SINR_averager.average(10.^(AWGN_SNR(:)/10),CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
            CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
            CQI(:,:,:) = CQI_temp*ones(LTE_params.Nrb,2,rank_loop);
        else  % single CQI value for every resource block
            for RB_ii = 1:LTE_params.Nrb
                if LTE_params.feedback.channel_averaging
                    rb_ii = RB_ii;
                else
                    rb_ii = (RB_ii-1)*LTE_params.Nsc+1:RB_ii*LTE_params.Nsc;
                end
                for slot_ii = 1:2
%                     if rank_i == 4
%                         AWGN_SNR = squeeze(SNR_g(i2,rb_ii,slot_ii,(i2-1)*2+1:min(i2*2,rank_i)));
%                     elseif rank_i == 3 && i2 == 2
%                         AWGN_SNR = squeeze(SNR_g(rank_i,rb_ii,slot_ii,i2:rank_i));
%                     else
                        AWGN_SNR = squeeze(SNR_g(rb_ii,slot_ii,1:rank_i));
%                     end
                    CQIs = (0:15);
                    SINReff = UE.SINR_averager.average(10.^(AWGN_SNR(:)/10),CQIs,[LTE_params.CQI_params(20).modulation_order,LTE_params.CQI_params(CQIs(2:end)).modulation_order]);
                    CQI_temp = LTE_common_CQI_mapping_table(LTE_params.CQI_mapping,SINReff,CQIs+1);
                    CQI(RB_ii,slot_ii,:) = CQI_temp*ones(1,1,rank_loop);
                end
            end
        end
%     end
else
    CQI = [];
end
