classdef MaxMinScheduler < network_elements.lteScheduler
% Scheduler that tries to maximize the minimum user throughput
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
       function obj = MaxMinScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets)
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
           N_UE = size(UE_output,2);
           N_RB = size(UE_output(1).CQI,1)*2;
 
           %% set pmi and ri values
           [UE_scheduling,c,user_ind] = obj.set_pmi_ri(UE_scheduling,N_UE,N_RB,UE_output);   
           count = 0;
%            c = zeros(N_RB,N_UE);
%            for u_= user_ind % calculate sum efficiencies for both codewords on every RB
%                count = count+1;
%                UE_scheduling(u_).UE_mapping = false(obj.RB_grid_size,2);
%                c(:,count) = sum(reshape([obj.CQI_params(UE_output(u_).CQI(:,:,1:UE_scheduling(u_).nCodewords)).efficiency],obj.RB_grid_size*2,UE_scheduling(u_).nCodewords),2);
%            end % based on this value the decision is carried out which UE gets an RB
%            c = c';
           C = zeros(N_UE,N_RB*N_UE);
           for rb = 1:N_RB
               for uu = 1:N_UE
                    C(uu,N_UE*(rb-1)+uu) = c(uu,rb);
               end
           end
           
           %% Max. min. scheduler
           RBs = obj.Max_min_scheduler(N_UE,N_RB,C);

           %% set cqi values
           UE_scheduling = obj.set_cqi(UE_scheduling,user_ind,RBs(1:end-1),N_UE,N_RB,UE_output);
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms); 
       end
       
       function RBs = Max_min_scheduler(obj,N_UE,N_RB,C)
           % core scheduling function (same in LL and SL)
           A = kron(eye(N_RB),ones(1,N_UE));
           A = [A,zeros(size(A,1),1)];
           constraints = [-C,-ones(size(C,1),1)];
           constraints = [constraints;A];
           RBs = linprog([zeros(N_UE*N_RB,1);1],constraints,[zeros(N_UE,1);ones(N_RB,1)],[],[],[zeros(N_RB*N_UE,1);-Inf],[ones(N_RB*N_UE,1);Inf],[],obj.linprog_options);
           for rb = 1:N_RB % randomize UE choice if RBs is noninteger
               toss = rand;
               RB_temp = RBs((rb-1)*N_UE+1:rb*N_UE);
               temp_ind = cumsum(RB_temp) > toss;
               RB_temp = zeros(size(RB_temp));
               RB_temp(find(temp_ind == 1,1)) = 1;
               RBs((rb-1)*N_UE+1:rb*N_UE) = RB_temp;
           end
%            RBs = IP([zeros(N_UE*N_RB,1);1],constraints,[zeros(N_UE,1);ones(N_RB,1)],[],[],[zeros(N_RB*N_UE,1);-Inf],[ones(N_RB*N_UE,1);Inf],1:24,10^-1);
       end
   end
end 
