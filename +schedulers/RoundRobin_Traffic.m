classdef RoundRobin_Traffic < network_elements.lteScheduler
% Round Robin scheduler that supports traffic models
% Martin Müller
% (c) 2011 by ITC
% www.nt.tuwien.ac.at

   properties
       SINR_averager
       CQI_mapping_data
       linprog_options
       alphabets
       av_throughput % exponentially weighted throughputs
       Tsub     % subframe duration
       alpha    % alpha parameter for alpha-fair utility functions
       J        % Jain's fairness index
       cqi_pdf            % empirical cqi pdf for each user - necessary to compute the right alpha value
       fairness           % desired fairness (for var. fair scheduler)
       efficiencies       % stores the spectral efficiencies corresponding to each cqi
       counter
       TP_exp
       last_set
       lambda
       rate_constraints
       lambda_store;
       MUMIMO
   end

   methods
       function obj = RoundRobin_Traffic(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets)
           % Fill in basic parameters (handled by the superclass constructor)
           obj = obj@network_elements.lteScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params);

           obj.static_scheduler = false;
           obj.SINR_averager = averager;
           % Get a vector of scheduling params (one for each UE)
           % initialized to the values that we want
           obj.UE_static_params = obj.get_initialized_UE_params(scheduler_params,CQI_params);
           obj.CQI_mapping_data = mapping_data;
           obj.alphabets = alphabets;
           obj.av_throughput = zeros(size(UEs_to_be_scheduled,2),1);
           obj.alpha = scheduler_params.alpha;
           obj.J = 1/size(UEs_to_be_scheduled,2);
           obj.cqi_pdf = zeros(size(UEs_to_be_scheduled,2),16);
           obj.fairness = scheduler_params.fairness;
           obj.efficiencies = repmat([obj.CQI_params(1).efficiency/2,obj.CQI_params(1:15).efficiency],size(UEs_to_be_scheduled,2),1);
           obj.counter = 0;
           obj.TP_exp = zeros(size(UEs_to_be_scheduled,2),1);
           obj.last_set = -10^6;
           obj.lambda = ones(size(UEs_to_be_scheduled,2),1);
           obj.lambda_store = cell(size(UEs_to_be_scheduled,2),1);
           obj.rate_constraints = scheduler_params.rate_constraints;
           obj.MUMIMO = scheduler_params.MUMIMO;
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           UE_scheduling = obj.UE_static_params;
           N_RB = size(UE_output(1).CQI,1)*2;
           N_UE = size(UE_output,2);
           for u_ = 1:N_UE
                obj.UEs(u_).traffic_model.check_TTI;
           end

           pdf_tmp = zeros(N_UE,16);
           pdf_const = obj.av_const;
           for uu = 1:N_UE
               CQI_temp = mod(UE_output(uu).CQI(:),20)+1;
               for rb = 1:N_RB
                    pdf_tmp(uu,CQI_temp(rb)) = pdf_tmp(uu,CQI_temp(rb))+1;
               end
               pdf_tmp(uu,:) = pdf_tmp(uu,:)/sum(pdf_tmp(uu,:));
           end
           obj.counter = obj.counter+1;
           if obj.counter ~= 1
               obj.cqi_pdf = pdf_tmp*1/pdf_const+obj.cqi_pdf*(1-1/pdf_const);
           else
               obj.cqi_pdf = pdf_tmp;
           end

%             if N_RB > 12
%                 error('The adjustable fairness scheduler is just designed for 1.4 MHz, please generalize me to arbitrary bandwidths')
%             end

           %% set pmi and ri values
           [UE_scheduling,c,user_ind] = obj.set_pmi_ri(UE_scheduling,N_UE,N_RB,UE_output);

           %% update average throughput
           for uu = 1:N_UE
               [obj.av_throughput(uu)] = obj.compute_av_throughput(UE_output(uu),obj.av_throughput(uu),uu,N_UE);
           end
           %% VF scheduler
           RBs = obj.Constrained_scheduler(N_UE,UE_output,c,user_ind,subframe_corr,SyncUsedElements,PBCHsyms);

           %% set cqi values
           UE_scheduling = obj.set_cqi(UE_scheduling,user_ind,RBs,N_UE,N_RB,UE_output);
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);

      end

       function RBs = Constrained_scheduler(obj,N_UE,UE_output,c,user_ind,subframe_corr,SyncUsedElements,PBCHsyms)

           overhead = SyncUsedElements+reshape(PBCHsyms,size(SyncUsedElements))+obj.overhead_ref; % overhead symbols;
           RBs = zeros(N_UE*2*obj.RB_grid_size,2); % resource block assignement
           utility = zeros(size(c)); % utility of the users
           symbols = kron(ones(1,size(c,3)),(obj.Ns_RB-overhead(:)));
           symbols = shiftdim(symbols,-1); % number of symbols of the different resource blocks (overhead subtracted)

           PMI_tmp = cat(3,UE_output(user_ind).PMI); % Overall PMI feedback
           RI_tmp = [UE_output(user_ind).RI]; % Overall RI feedback
           if isempty(RI_tmp)  % if there is no RI feedback --> we set the static values
               for uu = 1:N_UE
                    RI_tmp(uu) = obj.UE_static_params(uu).nLayers;
               end
           end
           if isempty(PMI_tmp)  % if there is no PMI feedback --> we use SUMIMO
               obj.MUMIMO = false;
           end
%            min_rank = min(RI_tmp);
%            max_rank = max(RI_tmp);

           if ~obj.MUMIMO && size(c,3) > 1 % in SUMIMO the efficiencies of both streams are simply added
               c = sum(c,3);
               c(:,:,2:end) = 0;
           end

           RB_UE = zeros(N_UE,2*obj.RB_grid_size+N_UE-1,2);
           bits_left = zeros(1,N_UE);
           N_UE_tmp = N_UE;
           UE_ind = 1:N_UE;
           user_ind_tmp = user_ind;
           RB_mat = eye(N_UE);
           RB_mat_tmp = RB_mat;
           state = 0;
           
%            user_ind_fb = user_ind;
%            user_ind_tr = user_ind;
%            for ii = 1:N_UE
%                if strcmp(obj.UEs(user_ind(ii)).traffic_model,'fullbuffer')
%                    user_ind_tr(ii) = 0;
%                else
%                    user_ind_fb(ii) = 0;
%                end
%            end
%            
%            user_ind_tr = user_ind(user_ind_tr~=0);
%            user_ind_fb = user_ind(user_ind_fb~=0);
           
%            for nn = 1:N_UE
%                 fprintf('\t\tuser %d', nn)
%                 obj.UEs(nn).traffic_model.bit_count
%            end
%            obj.UEs(4).traffic_model.bit_count
%            balken = '___________________________________________'

%            for uu = 1:N_UE
%                 RI_tmp(uu) = obj.UE_static_params(uu).nLayers;
%            end


           for rr = 1:2*obj.RB_grid_size                    % loop to assign resources                                              
               if state == 0 && ~isempty(user_ind_tmp)      % this is entered again, if all users with data have been scheduled once and there bit_count has been decreased
                    for ii = 1:length(user_ind_tmp)         % Check if data is available
                        if strcmp(obj.UEs(user_ind_tmp(ii)).traffic_model.type,'fullbuffer')
                            bits_left(user_ind_tmp(ii)) = 1;
                        else                        
                            bits_left(user_ind_tmp(ii)) = obj.UEs(user_ind_tmp(ii)).traffic_model.bit_count;
                        end
                            if ~bits_left(user_ind_tmp(ii))     %if no data for a user left, he will not be scheduled anymore
                                user_ind_tmp(ii) = 0;
                                N_UE_tmp = N_UE_tmp-1;                            
                            end
                    end

                   user_ind_tmp = user_ind_tmp(user_ind_tmp ~= 0);  %new user-vector without the ones without data

                   RB_mat_tmp = RB_mat(:,user_ind_tmp);     %columns interchanged in the random order of user_ind and without users without data

                   RB_UE(:,rr:(rr+N_UE_tmp-1),1) = RB_mat_tmp;
                   
                   state = mod(1,N_UE_tmp);     %if N_UE_tmp == 0 --> state = 0; so with only user, he will be checked for data during the next iteration-step; (if one uses 'state = 1' this is not the case)
               else
                   state = mod(state + 1, N_UE_tmp);    %if state = N_UE_tmp, users will be checked for data in the next iteration-step
               end
               
               if sum(bits_left) && ~strcmp(obj.UEs(UE_ind(UE_ind'.*RB_UE(:,rr,1)~=0)).traffic_model.type,'fullbuffer')
                   if  rr <= N_UE_tmp     %coarse decrease with crc-bits subtracted (only for non-fullbuffer)
                        bits_left(UE_ind(UE_ind'.*RB_UE(:,rr,1)~=0)) = obj.UEs(UE_ind(UE_ind'.*RB_UE(:,rr,1)~=0)).traffic_model.coarse_decrease(c(UE_ind(UE_ind'.*RB_UE(:,rr,1)~=0),rr,1)*(symbols(rr)-24));
                   elseif sum(bits_left)    %crc is subtracted only once
                       bits_left(UE_ind(UE_ind'.*RB_UE(:,rr,1)~=0)) = obj.UEs(UE_ind(UE_ind'.*RB_UE(:,rr,1)~=0)).traffic_model.coarse_decrease(c(UE_ind(UE_ind'.*RB_UE(:,rr,1)~=0),rr,1)*symbols(rr));
                   end
               end
               
%                  if sum(bits_left)
%                      bits_left(UE_ind(UE_ind'.*RB_UE(:,rr)~=0)) = obj.UEs(UE_ind(UE_ind'.*RB_UE(:,rr)~=0)).traffic_model.coarse_decrease(c(UE_ind(UE_ind'.*RB_UE(:,rr)~=0),rr,1)*symbols(rr));
%                  end
           end
          
%            RBs(:,1) = reshape(RB_UEs(:,:,1),[],1);
%            RBs(:,2) = reshape(RB_UEs(:,:,2),[],1);
%            
%            if ~sum(RBs(:,2)) % just one layer is used
%                RBs = RBs(:,1);
%            end
%            sum(RB_UE,2)

           true_RB_UE = RB_UE(user_ind,1:2*obj.RB_grid_size,1);    %putting the assigned_RBs in the right order, according to random order of user_ind

           RBs(:,1) = reshape(true_RB_UE,[],1);     %the size of true_RB_UE may be larger than the actual RB_grid; here it is truncated; if the number of RBs is not a multiple of the users with data, only the first few are scheduled; due to the random assignement this will change for every TTI
           
           RBs = RBs(:,1);
       end
   end
end