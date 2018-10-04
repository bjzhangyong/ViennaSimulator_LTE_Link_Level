classdef Kwan_scheduler < network_elements.lteScheduler
% Scheduler that tries to maximize the overall throughput based on
% "Multiuser Scheduling on the Downlink of an LTE Cellular System"
% R. Kwan, C. Leung, J. Zhang, Research Letters in Communications 2008
% Stefan Schwarz 
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

   properties
       SINR_averager
       CQI_mapping_data
       linprog_options
       alphabets
   end

   methods
       function obj = Kwan_scheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets)
           % Fill in basic parameters (handled by the superclass constructor)
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
           obj.linprog_options = optimset('LargeScale','off','Simplex','on','Display','off');
           obj.alphabets = alphabets;
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           UE_scheduling = obj.UE_static_params;
           N_RB = size(UE_output(1).CQI,1)*2;
           N_UE = size(UE_output,2);
           for uu = 1:N_UE
               UE_scheduling(uu).CQI_params = [];
               if ~isempty(UE_output(uu).PMI)   % check wheter PMI is fed back
                   if ~isempty(UE_output(uu).RI)    % check wheter RI is fed back
                        UE_scheduling(uu).nLayers = UE_output(uu).RI;
                        UE_scheduling(uu).nCodewords = min(2,UE_output(uu).RI);
                   end
               UE_scheduling(uu).PMI = UE_output(uu).PMI;
               end
           end
           % first stage: determine RB sets assigned to each UE
           % according to the CQI
           for u_=1:N_UE % calculate sum efficiencies for both codewords on every RB
               UE_scheduling(u_).UE_mapping = false(obj.RB_grid_size,2);
               eff_temp(:,u_) = sum(reshape([obj.CQI_params(UE_output(u_).CQI(:,:,1:UE_scheduling(u_).nCodewords)).efficiency],obj.RB_grid_size*2,UE_scheduling(u_).nCodewords),2);
           end % based on this value the decision is carried out which UE gets an RB
%            [~,Indis] = max(eff_temp,[],2);

           Indis = zeros(size(eff_temp,1),1);
           for r_ = 1:size(eff_temp,1)  % include randomization into the decision to make the scheduler fair if some UEs have equal CQIs
              maxi = max(eff_temp(r_,:));
              indis = find(eff_temp(r_,:) == maxi);
              ind = randi(length(indis),1);
              Indis(r_) = indis(ind);
           end
           for u_=1:N_UE
               UE_scheduling(u_).UE_mapping(Indis == u_) = true; 
           end
%            UE_output.CQI
%            UE_scheduling.UE_mapping
           % second stage: determine rate for each UE
           rate = zeros(N_UE,15);
           rb_num = zeros(N_UE,15);
           for uu = 1:N_UE
               for cqi_i = 1:15
                   temp_rate = 0;
                   rb_temp = 0;
                   for rb = 1:N_RB
                       if UE_scheduling(uu).UE_mapping(rb)
                           if UE_output(uu).CQI(rb) >= cqi_i && UE_output(uu).CQI(rb) ~= 20
                               temp_rate = temp_rate + obj.CQI_params(cqi_i).efficiency;
                               rb_temp = rb_temp+1;
                           elseif UE_output(uu).CQI(rb) == 20 && cqi_i == 1
                               temp_rate = temp_rate + obj.CQI_params(1).efficiency;
                               rb_temp = rb_temp+1;
                           end
                       end
                   end
                   rate(uu,cqi_i) = temp_rate;
                   rb_num(uu,cqi_i) = rb_temp;
               end
           end
           for uu = 1:N_UE
               [~,CQI_UE] = max(rate(uu,:));
               UE_output(uu).CQI(UE_output(uu).CQI == 20) = 1;
               UE_scheduling(uu).UE_mapping(UE_output(uu).CQI < CQI_UE) = 0;
               UE_scheduling(uu).cqi = CQI_UE;
               UE_scheduling(uu).assigned_RBs = squeeze(sum(sum(UE_scheduling(uu).UE_mapping,1),2));
               UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi,obj.CQI_params)];
           end
%            UE_scheduling.UE_mapping
%            UE_scheduling.cqi
%            error('jetzt')
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
       end
   end
end 
