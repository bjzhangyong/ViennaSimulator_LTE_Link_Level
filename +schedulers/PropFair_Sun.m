classdef PropFair_Sun < network_elements.lteScheduler
% OFDMA Proportional Fair Scheduler based on a greedy algorithm presented in
% "Reduced-Complexity Proportional Fair Scheduling for OFDMA Systems"
% Z. Sun, C. Yin, G. Yue, IEEE 2006 International Conference on Communications, Circuits and Systems Proceedings 
% Stefan Schwarz 
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

   properties
       SINR_averager
       CQI_mapping_data
       linprog_options
       alphabets
       av_throughput % exponentially weighted throughputs
       Tsub     % subframe duration
   end

   methods
       function obj = PropFair_Sun(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets,Tsub)
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
           obj.av_throughput = zeros(size(UEs_to_be_scheduled,2),1);
           obj.Tsub = Tsub;
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           UE_scheduling = obj.UE_static_params;
           N_RB = size(UE_output(1).CQI,1)*2;
           N_UE = size(UE_output,2);
           
           %% set pmi and ri values
           [UE_scheduling,c,user_ind] = obj.set_pmi_ri(UE_scheduling,N_UE,N_RB,UE_output);          
           c = c';
           
           %% update average throughput
           for uu = 1:N_UE
               obj.av_throughput(uu) = obj.compute_av_throughput(UE_output(uu),obj.av_throughput(uu),uu,N_UE);
           end
            
           %% PF scheduler
           RBs = obj.PF_scheduler(N_UE,N_RB,c,user_ind);
           
           %% set cqi values
           UE_scheduling = obj.set_cqi(UE_scheduling,user_ind,RBs,N_UE,N_RB,UE_output);
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);                     
       end
       
       function RBs = PF_scheduler(obj,N_UE,N_RB,c,user_ind)
           % core scheduling function (same in LL and SL)
           RB_set = ones(N_RB,1);
           RB_UEs = false(N_RB,N_UE);
           for rr = 1:N_RB
               res = find(RB_set);
               metric = zeros(length(res),N_UE);
               for r_ = 1:sum(RB_set)
                   for u_ = 1:N_UE
                       metric(res(r_),u_) = c(res(r_),u_)*12*7/((obj.av_const-1)*obj.av_throughput(user_ind(u_))+RB_UEs(:,u_).'*c(:,u_)*12*7);      % 12*7 equals the number of elements in a RB             
                   end
               end
               maxi = max(metric(:));
               indis = find(metric == maxi);
               ind = indis(randi(length(indis)));
               [temp_res,temp_ue] = ind2sub(size(metric),ind);
               RB_set(temp_res) = 0;
               RB_UEs(temp_res,temp_ue) = true;
           end
           RBs = vec(RB_UEs');
       end
   end
end 
