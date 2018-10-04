classdef MaxThroughputScheduler2 < network_elements.lteScheduler
% Scheduler that tries to maximize the overall throughput
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
       function obj = MaxThroughputScheduler2(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets)
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
           obj.linprog_options = optimset('LargeScale','off','Simplex','off','Display','off');
           obj.alphabets = alphabets;
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           UE_scheduling = obj.UE_static_params;
           c =[];
           N_UE = size(UE_output,2);
           N_RB = size(UE_output(1).CQI,1)*2;
%            permut = randperm(N_UE); % random permutation of the UEs (otherwise the linear programming resource block assignment is not fair)
           a = zeros(N_UE,1);
           b = zeros(N_UE,1);
           for uu = 1:N_UE
               UE_scheduling(uu).CQI_params = [];
               if ~isempty(UE_output(uu).PMI)   % check wheter PMI is fed back
                   if ~isempty(UE_output(uu).RI)    % check wheter RI is fed back
                        UE_scheduling(uu).nLayers = UE_output(uu).RI;
                        UE_scheduling(uu).nCodewords = min(2,UE_output(uu).RI);
                   end
               UE_scheduling(uu).PMI = UE_output(uu).PMI;
               end
               temp_UE(uu).CQI = UE_output(uu).CQI;
               temp_UE(uu).CQI(temp_UE(uu).CQI==20) = 0;
               CQI_bar = max(temp_UE(uu).CQI(:));
               if CQI_bar == 20
                   CQI_bar = 0;
               end
               CQI_bar = CQI_bar +1;
               a(uu) = obj.CQI_mapping_data.coeffs(1) * obj.SINR_averager.MI_data(1).k(CQI_bar);
               b(uu) = obj.CQI_mapping_data.coeffs(2) * obj.SINR_averager.MI_data(1).k(CQI_bar)+obj.SINR_averager.MI_data(1).d(CQI_bar);
%                UE_output(uu).CQI
           end
           for rb = 1:N_RB
               for uu = 1:N_UE
                    c = [c;a(uu)* temp_UE(uu).CQI(rb)+b(uu)];
               end
           end  
           
           % use best cqi to calculate a starting point for the linear
           % program --> speed increase
           for u_=1:obj.nUEs % calculate sum efficiencies for both codewords on every RB
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
           x0 = zeros(N_UE*N_RB,1);
           for u_ = 1:N_UE
               temp = zeros(size(Indis));
               temp(Indis == u_) = 1;
               x0(u_:N_UE:end) = temp;
           end
%            A = spdiags(ones(N_RB,N_UE),[0:N_UE-1],N_RB,N_UE+N_RB-1); 
           A = kron(eye(N_RB),ones(1,N_UE));
%            full(A)
           RBs = linprog(-c,A,ones(N_RB,1),[],[],zeros(N_RB*N_UE,1),ones(N_RB*N_UE,1),x0,obj.linprog_options);
           for uu = 1:N_UE
              temp_mapping = full(RBs(uu:N_UE:end));
              UE_scheduling(uu).UE_mapping = logical(reshape(temp_mapping,N_RB/2,2));
              UE_scheduling(uu).assigned_RBs = squeeze(sum(sum(obj.UE_static_params(uu).UE_mapping,1),2));
              UE_scheduling(uu).cqi = zeros(1,UE_scheduling(uu).nCodewords);
%               UE_scheduling(uu).UE_mapping
              if UE_scheduling(uu).assigned_RBs ~= 0
                  for i1 = 1:UE_scheduling(uu).nCodewords
                        CQI_all_temp = UE_output(uu).CQI(:,:,i1);
                        CQI_temp = CQI_all_temp(UE_scheduling(uu).UE_mapping);
                        SINRs = zeros(size(obj.CQI_mapping_data));
        %                             SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data(mod(CQI_temp,20)+1))/10),0:15);
                        SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data.table(mod(CQI_temp,20)+1)+0.9*obj.CQI_mapping_data.table(mod(CQI_temp,20)+2))/20),0:15,obj.alphabets); % this version is less conservative
                        SINRs = SINR_temp(:);
                        temp = zeros(size(SINRs));
                        temp(obj.CQI_mapping_data.table(1:16) <= SINRs) = 1;
                        temp_CQI = find(temp,1,'last')-1;
                        if temp_CQI
                            UE_scheduling(uu).cqi(i1) = temp_CQI;
                        else
                            UE_scheduling(uu).cqi(i1) = 20; % this is the rate 0 CQI
                        end
                        if logical(~UE_scheduling(uu).cqi(i1)) || (UE_scheduling(uu).cqi(i1) == 20)
                             UE_scheduling(uu).cqi(i1) = 1;
                        end
                        UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi(i1),obj.CQI_params)];
                  end 
              end
           end
%             UE_scheduling.UE_mapping
%             UE_scheduling.cqi
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
       end
   end
end 
