classdef VariableFairnessScheduler < network_elements.lteScheduler
% Scheduler that tries to maximize the throughput while guaranteeing a
% prespecified fairness level
% Stefan Schwarz 
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

   properties
       SINR_averager
       CQI_mapping_data
       linprog_options
       alphabets
       av_throughput
       fairness
   end

   methods
       function obj = VariableFairnessScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets)
           % Fill in basic parameters (handled by the superclass constructor)
            % Fill in basic parameters (handled by the superclass constructor)
           obj = obj@network_elements.lteScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params);
           
           % 0 CQI is not assigned (it means that conditions are out of
           % range to transmit, so no data will be assigned to 0-CQI RBs)
           %obj.assign_zero_CQI = scheduler_params.assign_zero_CQI;
           
           obj.static_scheduler = false;
           obj.SINR_averager = averager;
           obj.fairness = scheduler_params.fairness;
           % Get a vector of scheduling params (one for each UE)
           % initialized to the values that we want
           obj.UE_static_params = obj.get_initialized_UE_params(scheduler_params,CQI_params);
           obj.CQI_mapping_data = mapping_data;
           obj.linprog_options = optimset('LargeScale','off','Simplex','on','Display','off');
           obj.alphabets = alphabets;
           obj.av_throughput = zeros(length(UEs_to_be_scheduled),1);
           cvx_quiet(true); % surpress cvx output
           cvx_precision low;
%            cvx_solver sdpt3;
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           UE_scheduling = obj.UE_static_params;
           N_UE = size(UE_output,2);
           N_RB = size(UE_output(1).CQI,1)*2;
           
           %% set pmi and ri values
           [UE_scheduling,c,user_ind] = obj.set_pmi_ri(UE_scheduling,N_UE,N_RB,UE_output); 
           c = c*12*7; % convert from bits/cu to bits/RB
           C = zeros(N_UE,N_RB*N_UE);
           for rb = 1:N_RB
               for uu = 1:N_UE
                    C(uu,N_UE*(rb-1)+uu) = c(uu,rb);
               end
           end
           
           %% update average throughput
           for uu = 1:N_UE
               obj.av_throughput(uu) = obj.compute_av_throughput(UE_output(uu),obj.av_throughput(uu),uu,N_UE);
           end
 
           %% Variable fairness scheduler
           RBs = obj.var_fair_scheduler(N_UE,N_RB,c,C,user_ind);
           
           %% set cqi values
           UE_scheduling = obj.set_cqi(UE_scheduling,user_ind,RBs,N_UE,N_RB,UE_output);
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
                  end
       
       function RBs = var_fair_scheduler(obj,N_UE,N_RB,c,C,user_ind)
           % core scheduling function (same in LL and SL)
           F = kron(eye(N_RB),ones(1,N_UE));
           A = zeros(N_UE,N_RB*N_UE);
           for uu = 1:N_UE
                A(uu,uu:N_UE:end) = c(uu,:);
           end
           A = A*1/obj.av_const;
           b = (1-1/obj.av_const)*obj.av_throughput(user_ind);
           
           temp_fairness = obj.fairness;       
           c = vec(c); 
           d = sum(b);

           %% call to cvx to solve the SOCP
           cvx_begin
                variable RBs(N_UE*N_RB);
                maximize(c'*RBs);
                subject to
                    RBs <= 1;
                    RBs >= 0;
                    F*RBs == 1;    
                    norm(A*RBs+b) <= 1/sqrt(N_UE*temp_fairness)*(1/obj.av_const * c' * RBs + d);
           cvx_end
           if cvx_optval == -Inf
               disp('Prescribed fairness cannot be achieved - switching to max. min. fairness')
               A = kron(eye(N_RB),ones(1,N_UE));
               A = [A,zeros(size(A,1),1)];
               constraints = [-C,-ones(size(C,1),1)];
               constraints = [constraints;A];
               RBs = linprog([zeros(N_UE*N_RB,1);1],constraints,[zeros(N_UE,1);ones(N_RB,1)],[],[],[zeros(N_RB*N_UE,1);-Inf],[ones(N_RB*N_UE,1);Inf],[],obj.linprog_options);
               RBs = RBs(1:end-1);
           end
%            RBs = round(RBs);
           for rb = 1:N_RB % randomize UE choice if RBs is noninteger
               toss = rand;
               RB_temp = RBs((rb-1)*N_UE+1:rb*N_UE);
               temp_ind = cumsum(RB_temp) > toss;
               RB_temp = zeros(size(RB_temp));
               RB_temp(find(temp_ind == 1,1)) = 1;
               RBs((rb-1)*N_UE+1:rb*N_UE) = RB_temp;
           end
       end
       
   end
end 
