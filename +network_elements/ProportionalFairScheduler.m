classdef ProportionalFairScheduler < network_elements.lteScheduler
% Scheduler that dynamically adjusts PMI,RI and CQI and schedules users
% proportional to their current theoretically attainable rate
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

   properties
       SINR_averager
       CQI_mapping_data
   end

   methods
       function obj = ProportionalFairScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data)

           % Fill in basic parameters (handled by the superclass constructor)
           obj = obj@network_elements.lteScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params);
           
           % 0 CQI is not assigned (it means that conditions are out of
           % range to transmit, so no data will be assigned to 0-CQI RBs)
           %obj.assign_zero_CQI = scheduler_params.assign_zero_CQI;
           
           obj.static_scheduler = false;
           obj.SINR_averager = averager;
           % Get a vector of scheduling params (one for each UE)
           % initialized to the values that we want
           obj.UE_static_params = obj.get_initialized_UE_params(scheduler_params,CQI_params);
           obj.CQI_mapping_data = mapping_data;
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           % this function dynamically adjusts the transmission rank,
           % precoding matrix indicator and channel quality indicator
           % according to the feedback (if present); afterwards it schedules
           % users proportional to their theoretically attainable rate (as
           % the true one is not known)
           
%            UE_scheduling = obj.UE_static_params;
           UE_scheduling = obj.UE_static_params;
           for uu = 1:size(UE_output,2)
               UE_scheduling(uu).CQI_params = [];
               if ~isempty(UE_output(uu).PMI)   % check wheter PMI is fed back
                   if ~isempty(UE_output(uu).RI)    % check wheter RI is fed back
                        UE_scheduling(uu).nLayers = UE_output(uu).RI;
                        UE_scheduling(uu).nCodewords = min(2,UE_output(uu).RI);
                   end
               UE_scheduling(uu).PMI = UE_output(uu).PMI;
               end
%                if ~isempty(UE_output(uu).CQI)
%                    UE_scheduling(uu).cqi = squeeze(UE_output(uu).CQI);
% %                    if ~UE_scheduling(uu).cqi
% %                         UE_scheduling(uu).cqi = 1;
% %                    end
%                    for i1 = 1:UE_scheduling(uu).nCodewords
%                         UE_scheduling(uu).CQI_params(i1) = LTE_common_get_CQI_params(UE_scheduling(uu).cqi(i1),obj.CQI_params);
%                    end
%                end
           end

           rates = zeros(obj.nUEs,1);   % theoretical rates 
           for u_=1:obj.nUEs % calculate maximum rates, if different rates on different RBs were possible
           UE_scheduling(u_).UE_mapping = false(obj.RB_grid_size,2);
               for i_ = 1:UE_scheduling(u_).nCodewords
                   if ~isempty(UE_output(uu).CQI)
                       for r_ = 1:obj.RB_grid_size
                            for s_ = 1:2
                                rates(u_) = rates(u_) + obj.CQI_params(UE_output(u_).CQI(r_,s_,i_)).modulation_order*obj.Ns_RB*obj.CQI_params(UE_output(u_).CQI(r_,s_,i_)).coding_rate_x_1024 / 1024;
                            end
                       end
                   else
                       rates(u_) = rates(u_) + UE_scheduling(uu).CQI_params(i_).modulation_order*obj.Ns_RB*UE_scheduling(uu).CQI_params(i_).coding_rate_x_1024 / 1024 * 2*obj.RB_grid_size;
                       UE_output(u_).CQI(:,:,i_) = UE_scheduling(uu).CQI_params(i_).CQI * ones(obj.RB_grid_size,2);     % for later usage also if feedback is not available
                   end
               end 
               eff_temp(:,u_) = sum(reshape([obj.CQI_params(UE_output(u_).CQI(:,:,1:UE_scheduling(u_).nCodewords)).efficiency],obj.RB_grid_size*2,UE_scheduling(u_).nCodewords),2);
           end
           R = ones(obj.nUEs,1);   % real rates
           RB_set = true(obj.RB_grid_size*2,1);   % set of available resource blocks
           UE_set = ones(obj.nUEs,1);   % set of UEs that are evaluated
           CQI_set = zeros(obj.nUEs,2); % set of averaged CQI values
%            [~,max_u] = max(rates);
%            data = sum(reshape([obj.CQI_params(UE_output(max_u).CQI).modulation_order]*obj.Ns_RB.*[obj.CQI_params(UE_output(max_u).CQI).coding_rate_x_1024] / 1024,obj.RB_grid_size,2,UE_scheduling(max_u).nCodewords),3);
%            [R(max_u),max_RB] = max(data(:));
%            RB_set(max_RB) = 0;
           eff_temp = eff_temp + eps;
           rates = rates + eps;
           [mini,min_u] = min(R.*UE_set./rates);
           R = zeros(size(R));
           for k_ = 1:sum(sum(RB_set))
               if mini == Inf
                   break
               end
%                data = sum(reshape([obj.CQI_params(UE_output(min_u).CQI(:,:,1:UE_scheduling(min_u).nCodewords)).modulation_order]*obj.Ns_RB.*[obj.CQI_params(UE_output(min_u).CQI(:,:,1:UE_scheduling(min_u).nCodewords)).coding_rate_x_1024] / 1024,obj.RB_grid_size*2,UE_scheduling(min_u).nCodewords),2);   
%                [R_temp,max_RB] = max(data.*RB_set);
%                R(min_u) = R(min_u) + R_temp;
%                RB_set(max_RB) = false;
%                UE_scheduling(min_u).UE_mapping(max_RB) = true;
                CQI_indis = [];
                CQI_old = CQI_set(min_u,:);
                while 1
                   [temp_rate,CQI_ind,nomore,CQIs] = obj.get_rate(eff_temp(:,min_u),RB_set,UE_scheduling(min_u),UE_output(min_u));
                   if temp_rate <= R(min_u)
                       if nomore
                            UE_set(min_u) = Inf;
                            UE_scheduling(min_u).UE_mapping(CQI_indis) = false; 
                            RB_set(CQI_indis) = true;
                            CQI_set(min_u,:) = CQI_old;
                            break;
                       else
                           CQI_indis = [CQI_indis,CQI_ind];
                           UE_scheduling(min_u).UE_mapping(CQI_ind) = true;
                           RB_set(CQI_ind) = false;
                           CQI_set(min_u,:) = CQIs;
                       end
                   else
                       R(min_u) = temp_rate;
                       RB_set(CQI_ind) = false;
                       UE_scheduling(min_u).UE_mapping(CQI_ind) = true; 
                       CQI_set(min_u,1:UE_scheduling(min_u).nCodewords) = CQIs';
                       break;
                   end
                end
               [mini,min_u] = min(R.*UE_set./rates);
           end
           for u_ = 1:obj.nUEs
                UE_scheduling(u_).assigned_RBs = squeeze(sum(sum(obj.UE_static_params(u_).UE_mapping,1),2));
                if UE_scheduling(u_).assigned_RBs
                    for i_ = 1:UE_scheduling(u_).nCodewords
                        UE_scheduling(u_).CQI_params = [UE_scheduling(u_).CQI_params,LTE_common_get_CQI_params(CQI_set(u_,i_),obj.CQI_params)];
                        UE_scheduling(u_).cqi(i_) = CQI_set(u_,i_);
                    end
                end
           end
           
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
       end
       
       function [temp_rate,CQI_ind,nomore,CQIs] = get_rate(obj,eff_temp,RB_set,UE_scheduling,UE_output)
           if sum(RB_set) == 0
               nomore = true;
               temp_rate = 0;
               CQI_ind = 0;
               CQIs = 0;
           else
               CQIs = zeros(1,UE_scheduling.nCodewords);
               [~,CQI_ind] = max(eff_temp.*RB_set);
               UE_scheduling.UE_mapping(CQI_ind) = true;
               for c_ = 1:UE_scheduling.nCodewords
                    mapping = false([size(UE_scheduling.UE_mapping),UE_scheduling.nCodewords]);
                    mapping(:,:,c_) = UE_scheduling.UE_mapping;
                    CQI_temp = UE_output.CQI(mapping);
                    SINRs = zeros(size(obj.CQI_mapping_data));
%                     SINR_temp = obj.SINR_averager.average(10.^(obj.CQI_mapping_data(mod(CQI_temp,20)+1)/10),0:15);
                    SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data(mod(CQI_temp,20)+1)+obj.CQI_mapping_data(mod(CQI_temp,20)+2))/20),0:15); % this version is less conservative
                    SINRs = SINR_temp(:);
                    temp = zeros(size(SINRs));
                    temp(obj.CQI_mapping_data(1:16) <= SINRs) = 1;
                    temp_CQI = find(temp,1,'last')-1;
                    if temp_CQI
                        CQIs(1,c_) = temp_CQI;
                    else
                        CQIs(1,c_) = 20; % this is the rate 0 CQI
                    end
%                     CQIs(1,c_) = find(temp,1,'last');
               end
               temp_rate = sum([obj.CQI_params(CQIs).modulation_order]*sum(sum(UE_scheduling.UE_mapping)).*[obj.CQI_params(CQIs).coding_rate_x_1024] / 1024)+eps*sum(sum(UE_scheduling.UE_mapping)); % the eps term is necessary to set a RB if the SNR is very weak
               CQIs(CQIs == 20) = 1;  % CQI 0 (=20) is not possible in coding
               nomore = false;
           end
       end
           
   end
end 
