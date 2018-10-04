classdef Weighted_sumrate < network_elements.lteScheduler
% Weighted sum rate maximization scheduler
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
       weights  % weights for the weighted sum rate maximization
       cqi_pmf            % empirical cqi pmf for each user - necessary to compute the appropriate weights
       cqi_pmf_lin
       counter = 0;
       efficiencies
       rate_pmf
       exp_rate
       rate_constraints
       rate_store = [];
   end

   methods
       function obj = Weighted_sumrate(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets,weights)
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
           obj.weights = weights;
%            obj.cqi_pmf = zeros(size(UEs_to_be_scheduled,2),obj.RB_grid_size*2,16);
           obj.cqi_pmf = zeros(size(UEs_to_be_scheduled,2),16);
           obj.counter = 0;
           obj.efficiencies = repmat([obj.CQI_params(20).efficiency,obj.CQI_params(1:15).efficiency],size(UEs_to_be_scheduled,2),1)*(1-1/2.5)+repmat([obj.CQI_params(1:15).efficiency,obj.CQI_params(15).efficiency],size(UEs_to_be_scheduled,2),1)*1/2.5;
%            obj.rate_pmf = zeros(size(UEs_to_be_scheduled,2),16);
           obj.rate_pmf = ones(size(UEs_to_be_scheduled,2),1);
           obj.cqi_pmf_lin = zeros(size(UEs_to_be_scheduled,2),16);
           obj.rate_constraints = scheduler_params.rate_constraints;
           obj.exp_rate = zeros(size(UEs_to_be_scheduled,2),1);
       end

      function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           UE_scheduling = obj.UE_static_params;
           N_UE = size(UE_output,2);
           N_RB = size(UE_output(1).CQI,1)*2;
           
           pmf_tmp = zeros(size(obj.cqi_pmf));
           pmf_const = obj.av_const;
           obj.cqi_pmf_lin = obj.cqi_pmf_lin*obj.counter*N_RB;
           for uu = 1:N_UE
               CQI_temp = mod(UE_output(uu).CQI(:),20)+1;
               for rb = 1:N_RB
                    pmf_tmp(uu,CQI_temp(rb)) = pmf_tmp(uu,CQI_temp(rb))+1;
%                     obj.cqi_pmf_store(uu,CQI_temp(rb),end) = obj.cqi_pmf_store(uu,CQI_temp(rb),end)+1;
                    obj.cqi_pmf_lin(uu,CQI_temp(rb)) = obj.cqi_pmf_lin(uu,CQI_temp(rb))+1;
%                      pmf_tmp(uu,rb,:) = pmf_tmp(uu,rb,:)./sum(pmf_tmp(uu,rb,:));
               end
               pmf_tmp(uu,:) = pmf_tmp(uu,:)/sum(pmf_tmp(uu,:));
           end
           obj.counter = obj.counter+1;
           obj.cqi_pmf_lin = obj.cqi_pmf_lin/(obj.counter*N_RB);
           if obj.counter ~= 1
               obj.cqi_pmf = pmf_tmp*1/pmf_const+obj.cqi_pmf*(1-1/pmf_const);
           else
               obj.cqi_pmf = pmf_tmp;
           end  
           obj.cqi_pmf = obj.cqi_pmf_lin;
           
           %% set pmi and ri values
           [UE_scheduling,c,user_ind] = obj.set_pmi_ri(UE_scheduling,N_UE,N_RB,UE_output); 
           for uu = 1:N_UE
                if ~isempty(find(UE_output(uu).CQI(:) == 20,1))
                    obj.efficiencies(uu,1) = min(c(uu,:));
                end
           end
           for i = 1:N_UE
               c(i,:) = c(i,:)*obj.weights(user_ind(i));
           end
%            if ~mod(obj.counter,100)
% 
%             obj.weight_adaptation(N_UE);
% %             obj.expected_rate(N_UE);
%            end
           c = vec(c);
           
           %% Max sum rate scheduler
           RBs = obj.WSR_scheduler(N_UE,N_RB,c);
           
           %% set cqi values
           UE_scheduling = obj.set_cqi(UE_scheduling,user_ind,RBs,N_UE,N_RB,UE_output);
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
       end
       
       function RBs = WSR_scheduler(obj,N_UE,N_RB,c)
           % core scheduling function (same in LL and SL)
           A = kron(eye(N_RB),ones(1,N_UE));
           RBs = linprog(-c,A,ones(N_RB,1),[],[],zeros(N_RB*N_UE,1),ones(N_RB*N_UE,1),[],obj.linprog_options);  
           % call to cvx to solve the problem
%            cvx_begin
%                 variable RBs(N_UE*N_RB);
%                 maximize(c'*RBs);
%                 subject to
%                     RBs <= 1;
%                     RBs >= 0;
%                     A*RBs == 1;                
%            cvx_end
       end
       
%        function expected_rate(obj,N_UE)
%            c = zeros(N_UE,16);
%            Exp_TP = zeros(N_UE,16);
%            for uu = 1:N_UE
%                c(uu,:) = obj.efficiencies(uu,:)*obj.weights(uu);
%            end
%            for uu = 1:N_UE
% %                for rb = 1:obj.RB_grid_size*2
%                    for cqi_i = 1:16
%                        if obj.cqi_pmf(uu,cqi_i) ~= 0
%                            P_smaller = ones(N_UE-1,1);
%                            P_equal = zeros(N_UE-1,1);
%                            count1 = 0;
%                            equal_counter = 0;
%                            equal_ind = [];
%                            for u2 = 1:N_UE
%                                if u2 ~= uu
%                                    count1 = count1+1;
%                                    indis1 = c(u2,1:16) < c(uu,cqi_i);
%                                    indis2 = c(u2,1:16) == c(uu,cqi_i);
%                                    P_smaller(count1) = sum(obj.cqi_pmf(u2,logical(indis1)));
%                                    if sum(indis2) ~= 0
%                                        P_equal(count1) = obj.cqi_pmf(u2,logical(indis2));
%                                    end
%                                    if P_equal(count1)
%                                         equal_counter = equal_counter+1;
%                                         equal_ind = [equal_ind,count1];
%                                    end                           
%                                end
%                            end
%                            P_temp = prod(P_smaller);
%                            for p1 = 1:equal_counter
%                                temp_ind = nchoosek(equal_ind,p1);
%                                for p2 = 1:size(temp_ind,1)
%                                    temp_ind_all = true(size(P_smaller));
%                                    temp_ind_all(temp_ind(p2,:)) = false;
%                                    P_temp = P_temp + 1/(size(temp_ind,2)+1)*prod(P_equal(temp_ind(p2,:)))*prod(P_smaller(temp_ind_all));
%                                end
%                            end
%                            Exp_TP(uu,cqi_i) = obj.cqi_pmf(uu,cqi_i)*P_temp;
%                        else
%                            Exp_TP(uu,cqi_i) = 0;
%                        end
%                    end
% %                end
%            end
% %            obj.rate_pmf = sum(Exp_TP,3)/rb;
%            obj.rate_pmf = Exp_TP;
% %            obj.rate_pmf(:,1) = obj.rate_pmf(:,1)+(1-sum(obj.rate_pmf,2));
%            temp_rate1 = 8*round(1/8*obj.rate_pmf.*obj.efficiencies*(obj.Ns_RB-obj.overhead_ref)*12)-24;
%            temp_rate2 = 8*round(1/8*obj.rate_pmf.*obj.efficiencies*(obj.Ns_RB-obj.overhead_sync)*12)-24;
%            temp_rate1(temp_rate1 < 0) = 0;
%            temp_rate2(temp_rate2 < 0) = 0;
%            temp_rate = temp_rate1*4/5+temp_rate2*1/5;
%            obj.exp_rate = sum(temp_rate,2);
% %            exp(-2*(1000-obj.exp_rate(1))^2/((obj.efficiencies(1,16)*(obj.Ns_RB-obj.overhead_ref)*12)^2))-1
% %            exp(-2*(1000-obj.exp_rate(2))^2/((obj.efficiencies(1,16)*(obj.Ns_RB-obj.overhead_ref)*12)^2))-1
%        end

       function rate_back = expected_rate(obj,N_UE)
%            c = zeros(N_UE,16);
%            Exp_TP = zeros(N_UE,16);
%            for uu = 1:N_UE
%                c(uu,:) = obj.efficiencies(uu,:)*obj.weights(uu);
%            end
           rate_back = [];
           c = zeros(N_UE,12);
           for i = 1:500
               for uu = 1:N_UE
                   ind_max = find(obj.cqi_pmf(uu,:) ~= 0,1,'last');
                   ind_min = find(obj.cqi_pmf(uu,:) ~= 0,1,'first');
                   int_tmp = randi([ind_min,ind_max],1);
                   range = min(3,int_tmp-ind_min);
                   c(uu,:) = obj.efficiencies(uu,randsrc(1,12,[int_tmp-range:int_tmp;obj.cqi_pmf(uu,int_tmp-range:int_tmp)/sum(obj.cqi_pmf(uu,int_tmp-range:int_tmp))]))*obj.weights(uu);
               end
%            for uu = 1:N_UE
% %                for rb = 1:obj.RB_grid_size*2
%                    for cqi_i = 1:16
%                        if obj.cqi_pmf(uu,cqi_i) ~= 0
%                            P_smaller = ones(N_UE-1,1);
%                            P_equal = zeros(N_UE-1,1);
%                            count1 = 0;
%                            equal_counter = 0;
%                            equal_ind = [];
%                            for u2 = 1:N_UE
%                                if u2 ~= uu
%                                    count1 = count1+1;
%                                    indis1 = c(u2,1:16) < c(uu,cqi_i);
%                                    indis2 = c(u2,1:16) == c(uu,cqi_i);
%                                    P_smaller(count1) = sum(obj.cqi_pmf(u2,logical(indis1)));
%                                    if sum(indis2) ~= 0
%                                        P_equal(count1) = obj.cqi_pmf(u2,logical(indis2));
%                                    end
%                                    if P_equal(count1)
%                                         equal_counter = equal_counter+1;
%                                         equal_ind = [equal_ind,count1];
%                                    end                           
%                                end
%                            end
%                            P_temp = prod(P_smaller);
%                            for p1 = 1:equal_counter
%                                temp_ind = nchoosek(equal_ind,p1);
%                                for p2 = 1:size(temp_ind,1)
%                                    temp_ind_all = true(size(P_smaller));
%                                    temp_ind_all(temp_ind(p2,:)) = false;
%                                    P_temp = P_temp + 1/(size(temp_ind,2)+1)*prod(P_equal(temp_ind(p2,:)))*prod(P_smaller(temp_ind_all));
%                                end
%                            end
%                            Exp_TP(uu,cqi_i) = obj.cqi_pmf(uu,cqi_i)*P_temp;
%                        else
%                            Exp_TP(uu,cqi_i) = 0;
%                        end
%                    end
% %                end
%            end
% %            obj.rate_pmf = sum(Exp_TP,3)/rb;
%            obj.rate_pmf = Exp_TP;
% %            obj.rate_pmf(:,1) = obj.rate_pmf(:,1)+(1-sum(obj.rate_pmf,2));
%            temp_rate1 = 8*round(1/8*obj.rate_pmf.*obj.efficiencies*(obj.Ns_RB-obj.overhead_ref)*12)-24;
%            temp_rate2 = 8*round(1/8*obj.rate_pmf.*obj.efficiencies*(obj.Ns_RB-obj.overhead_sync)*12)-24;
%            temp_rate1(temp_rate1 < 0) = 0;
%            temp_rate2(temp_rate2 < 0) = 0;
%            temp_rate = temp_rate1*4/5+temp_rate2*1/5;
           [~,ind] = max(c,[],1);
           tmp_rate = zeros(N_UE,12);
           for uu = 1:N_UE
               tmp_rate(uu,ind == uu) = c(uu,ind == uu)/obj.weights(uu);
           end
           temp_rate1 = max(zeros(N_UE,1),8*round(1/8*mean(tmp_rate,2)*(obj.Ns_RB-obj.overhead_ref)*12)-24);
           temp_rate2 = max(zeros(N_UE,1),8*round(1/8*mean(tmp_rate,2)*(obj.Ns_RB-obj.overhead_sync-obj.overhead_ref)*12)-24);
           temp_rate = temp_rate1*4/5+temp_rate2*1/5;
           rate_back = [rate_back,sum(temp_rate,2)];
%            obj.exp_rate = (1-1/obj.av_const)*obj.exp_rate+1/obj.av_const * sum(temp_rate,2);
%            exp(-2*(1000-obj.exp_rate(1))^2/((obj.efficiencies(1,16)*(obj.Ns_RB-obj.overhead_ref)*12)^2))-1
%            exp(-2*(1000-obj.exp_rate(2))^2/((obj.efficiencies(1,16)*(obj.Ns_RB-obj.overhead_ref)*12)^2))-1
           end
           rate_back = mean(rate_back,2);
       end

       
       function weight_adaptation(obj,N_UE)
           rate_back = obj.expected_rate(N_UE);
           obj.exp_rate = rate_back;
           step = 10^-4;
           for uu = 1:length(obj.rate_constraints)
               if isinf(obj.rate_constraints(uu))
                   continue
               end
               count = 0;
               if obj.exp_rate(uu)-obj.rate_constraints(uu) < 0 
                   while abs(obj.exp_rate(uu)-obj.rate_constraints(uu)) > 10 && count < 0
                       count = count+1;
%                        if obj.weights(uu) == 1
%                            disp(['rate constraint number ' num2str(uu) ' not feasible']);
%                            break;
%                        end
                       obj.weights(uu) = max(0,min(1,obj.weights(uu)+ step*(obj.rate_constraints(uu)-obj.exp_rate(uu))));
                       ind = true(size(obj.weights));
                       ind(uu) = false;
                       if obj.weights(uu) ~= 1
                            obj.weights(ind) = obj.weights(ind)/sum(obj.weights(ind))*(1-obj.weights(uu));
                       else 
                           obj.weights(ind) = 0;
                       end
                       rate_tmp = obj.expected_rate(N_UE);
%                        obj.exp_rate = (1-1/obj.av_const)*obj.exp_rate+1/obj.av_const * rate_tmp;
                       obj.exp_rate = rate_tmp;
%                        obj.rate_store = [obj.rate_store,obj.exp_rate];
                   end
               end
           end
       end
       
%        function Wang_iterations(N_UE)
%            obj.cqi_pmf
%        end
   end
end 
