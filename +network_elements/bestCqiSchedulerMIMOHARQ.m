classdef bestCqiSchedulerMIMOHARQ < network_elements.lteScheduler
% Scheduler that dynamically adjusts PMI,RI and CQI and schedules users
% on their RB with the best CQI. This class supersedes the old best CQI
% scheduler. Patch that adds HARQ support for the single user case. DO NOT
% USE THIS SCHEDULER FOR MORE THAN ONE USER. The CQI assignment will then
% be wrong. Used for the validation of the HARQ modeling.
% Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

   properties
       SINR_averager
       CQI_mapping_data
       alphabets
   end

   methods
       function obj = bestCqiSchedulerMIMOHARQ(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,varargin)
       
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
           if ~isempty(varargin)
               obj.alphabets = varargin{1};
           else
               obj.alphabets = [];
           end
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           % this function dynamically adjusts the transmission rank,
           % precoding matrix indicator and channel quality indicator
           % according to the feedback (if present); afterwards it schedules
           % users proportional to their theoretically attainable rate (as
           % the true one is not known)

           % set PMI and RI values
           N_UE = size(UE_output,2);
           N_RB = size(UE_output(1).CQI,1)*2;
           UE_scheduling = obj.UE_static_params;
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
           
           %% scheduler according to the CQI
           eff_temp = zeros(N_RB,N_UE);
           for u_=1:N_UE % calculate sum efficiencies for both codewords on every RB
               UE_scheduling(u_).UE_mapping = false(obj.RB_grid_size,2);
               eff_temp(:,u_) = sum(reshape([obj.CQI_params(UE_output(u_).CQI(:,:,1:UE_scheduling(u_).nCodewords)).efficiency],obj.RB_grid_size*2,UE_scheduling(u_).nCodewords),2);
           end % based on this value the decision is carried out which UE gets an RB

           RBs = false(N_UE*N_RB,1);
           for r_ = 1:N_RB  % include randomization into the decision to make the scheduler fair if some UEs have equal CQIs
              maxi = max(eff_temp(r_,:));
              indis = find(eff_temp(r_,:) == maxi);
              ind = randi(length(indis),1);
              ind = indis(ind);
              RBs(N_UE*(r_-1)+1+ind-1) = true;
%               Indis(r_) = indis(ind);
           end
           
           %%
           UE_scheduling = obj.set_cqi(UE_scheduling,1:N_UE,RBs(1:end),N_UE,N_RB,UE_output);
           
           % HARQ handling. Only for SUSISO case. This is a small patch to
           % the Berst CQI MIMO scheduler to enable validaiton simulations
           % of the HARQ model.
           if obj.UE_specific(1).current_HARQ_process(1).rv_idx > 0
               UE_scheduling(1).cqi = obj.UE_specific(1).current_HARQ_process(1).cqi;
           end
               
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms); 
       end
       
   end
end 
