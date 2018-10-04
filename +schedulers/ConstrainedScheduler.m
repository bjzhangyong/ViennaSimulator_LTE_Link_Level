classdef ConstrainedScheduler < network_elements.lteScheduler
% Variable fair scheduler based on the proportional fair scheduler presented in 
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
       function obj = ConstrainedScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,alphabets)
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
           obj.lambda = zeros(size(UEs_to_be_scheduled,2),1);
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
           %% store the cqi pdf
% %            obj.cqi_pdf_lin = obj.cqi_pdf_lin*obj.counter*N_RB;
%            obj.cqi_pdf_store = circshift(obj.cqi_pdf_store,[0,0,-1]);
%            obj.cqi_pdf_store(:,:,end) = 0; 
           pdf_tmp = zeros(N_UE,16);
           pdf_const = obj.av_const;
           for uu = 1:N_UE
               CQI_temp = mod(UE_output(uu).CQI(:),20)+1;
               for rb = 1:N_RB
                    pdf_tmp(uu,CQI_temp(rb)) = pdf_tmp(uu,CQI_temp(rb))+1;
%                     obj.cqi_pdf_store(uu,CQI_temp(rb),end) = obj.cqi_pdf_store(uu,CQI_temp(rb),end)+1;
% %                     obj.cqi_pdf_lin(uu,CQI_temp(rb)) = obj.cqi_pdf_lin(uu,CQI_temp(rb))+1;
               end
               pdf_tmp(uu,:) = pdf_tmp(uu,:)/sum(pdf_tmp(uu,:));
           end
           obj.counter = obj.counter+1;
% %            obj.cqi_pdf_lin = obj.cqi_pdf_lin/(obj.counter*N_RB);
           if obj.counter ~= 1
               obj.cqi_pdf = pdf_tmp*1/pdf_const+obj.cqi_pdf*(1-1/pdf_const);
           else
               obj.cqi_pdf = pdf_tmp;
           end
%            obj.cqi_pdf_lin = sum(obj.cqi_pdf_store(:,:,max(1,end-obj.counter):end),3);
%            obj.cqi_pdf_lin = obj.cqi_pdf_lin/(min(obj.counter,obj.av_const)*N_RB);

            if N_RB > 12
                error('The adjustable fairness scheduler is just designed for 1.4 MHz, please generalize me to arbitrary bandwidths')
            end
               
           %% set pmi and ri values
           [UE_scheduling,c,user_ind] = obj.set_pmi_ri(UE_scheduling,N_UE,N_RB,UE_output);          
%            for uu = 1:N_UE
%                 if ~isempty(find(UE_output(uu).CQI(:) == 20,1))
%                     obj.efficiencies(uu,1) = min(c(uu,:));
%                 end
%            end
%            c = zeros(N_RB,N_UE);
%            count = 0;
%            for u_= user_ind % calculate sum efficiencies for both codewords on every RB
%                count = count+1;
%                UE_scheduling(u_).UE_mapping = false(obj.RB_grid_size,2);
%                c(:,count) = sum(reshape([obj.CQI_params(UE_output(u_).CQI(:,:,1:UE_scheduling(u_).nCodewords)).efficiency],obj.RB_grid_size*2,UE_scheduling(u_).nCodewords),2);
%            end % based on this value the decision is carried out which UE gets an RB

           %% update average throughput
           for uu = 1:N_UE
%                [obj.av_throughput(uu),obj.av_throughput_lin(uu,:)] = obj.compute_av_throughput(UE_output(uu),obj.av_throughput(uu),uu,N_UE,obj.av_throughput_lin(uu,:),obj.counter);
               [obj.av_throughput(uu)] = obj.compute_av_throughput(UE_output(uu),obj.av_throughput(uu),uu,N_UE);
           end
           %% update Lagrange multipliers
           obj.lambda_update(N_UE,UE_output);

           %% find compatible users for MUMIMO
%            obj.Compatible_users(N_UE,UE_output,c,user_ind);

           %% VF scheduler
           RBs = obj.Constrained_scheduler(N_UE,UE_output,c,user_ind,subframe_corr,SyncUsedElements,PBCHsyms);
           
           %% set cqi values
           UE_scheduling = obj.set_cqi(UE_scheduling,user_ind,RBs,N_UE,N_RB,UE_output);
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms); 
           
          end
       
%        function RBs = Constrained_scheduler(obj,N_UE,N_RB,c,user_ind,subframe_corr)
%            % core scheduling function (same in LL and SL)
%            RB_set = true(N_RB,1);
%            RB_UEs = false(N_RB,N_UE);
%            alpha_tmp = obj.alpha(end);
%            if ~mod(subframe_corr-1,5)
%                overhead = obj.overhead_ref+obj.overhead_sync;
%            else
%                overhead = obj.overhead_ref;
%            end
%            for rr = 1:N_RB
%                res = find(RB_set);
%                metric = ones(N_RB,N_UE)*-Inf;
%                for r_ = 1:sum(RB_set)
%                    for u_ = 1:N_UE
%                        if ~strcmp(obj.UEs(user_ind(u_)).traffic_model.type,'fullbuffer') 
%                             bitnr = obj.UEs(user_ind(u_)).traffic_model.bit_count;
%                             if  bitnr <= 0
%                                 metric(res(r_),u_) = -10^6;
%                             else
%                                 metric(res(r_),u_) = min(c(res(r_),u_)*obj.Ns_RB,bitnr)*obj.lambda(user_ind(u_));
%                             end
%                        else   
%                             metric(res(r_),u_) = c(res(r_),u_)*obj.Ns_RB*(max(obj.av_throughput(user_ind(u_)),eps)^alpha_tmp+obj.lambda(user_ind(u_)));  
%                        end
%                    end
%                end
%                maxi = max(metric(:));
%                indis = find(metric == maxi);
%                ind = indis(randi(length(indis)));
%                [temp_res,temp_ue] = ind2sub(size(metric),ind);
%                RB_set(temp_res) = false;
%                if maxi > -10^6
%                    RB_UEs(temp_res,temp_ue) = true;
%                else
%                    RB_UEs(temp_res,temp_ue) = false;
%                end
%                if ~strcmp(obj.UEs(user_ind(temp_ue)).traffic_model.type,'fullbuffer')
%                    if isinf(obj.UEs(user_ind(temp_ue)).traffic_model.bit_count)
%                        obj.UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_res,temp_ue)*(obj.Ns_RB-overhead),1);
%                    else
%                        obj.UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_res,temp_ue)*(obj.Ns_RB-overhead));
%                    end
%                end
%            end
%            RB_UEs = RB_UEs';
%            RBs = RB_UEs(:);
%        end  
       
       function lambda_update(obj,N_UE,UE_output)
%           for uu = 1:N_UE
%                TP = 0;
%                if ~isempty(UE_output(uu).rx_data_bits)
%                    for cc = 1:length(UE_output(uu).rx_data_bits)
%                          TP = TP + UE_output(uu).ACK(cc)*length(UE_output(uu).rx_data_bits{cc});
%                    end
%                    if ~isinf(obj.rate_constraints(uu))
%                         obj.lambda(uu) = max(0,obj.lambda(uu)-1/obj.av_const*(TP-obj.rate_constraints(uu)));
%                         obj.lambda_store{uu} = [obj.lambda_store{uu},obj.lambda(uu)];
%                    end
%                    if strcmp(obj.UEs(uu).traffic_model.type,'voip')
%                        obj.lambda(uu) = max(0,obj.lambda(uu)-1/obj.av_const*(obj.UEs(uu).traffic_model.delay_constraint-obj.UEs(uu).traffic_model.get_buffer_length/obj.UEs(uu).traffic_model.get_arrival_rate));
%                        obj.lambda_store{uu} = [obj.lambda_store{uu},obj.lambda(uu)];
%                    end
%                end
%            end
       end
       
       function RBs = Constrained_scheduler(obj,N_UE,UE_output,c,user_ind,subframe_corr,SyncUsedElements,PBCHsyms)
           overhead = SyncUsedElements+reshape(PBCHsyms,size(SyncUsedElements))+obj.overhead_ref; % overhead symbols; 
           RBs = zeros(N_UE*2*obj.RB_grid_size,2); % resource block assignement
           alpha_tmp = obj.alpha(end);  % fairness parameter
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
           min_rank = min(RI_tmp);
           max_rank = max(RI_tmp);
           
%            if ~obj.MUMIMO && size(c,3) > 1 % in SUMIMO the efficiencies of both streams are simply added
%                c = c(:,:,1)+c(:,:,2);
%                c(:,:,2) = 0;
%            end
           if ~obj.MUMIMO && size(c,3) > 1 % in SUMIMO the efficiencies of both streams are simply added
               c = sum(c,3);
               c(:,:,2:end) = 0;
           end
           for uu = 1:N_UE % compute user utility for all resources
               symbols_tmp = symbols;
               size(obj.av_throughput)
%                if RI_tmp(uu) == 3
%                    symbols_tmp(:,:,1) = symbols_tmp(:,:,1)*2; % rank 3: two layers for codeword 1
%                elseif RI_tmp(uu) == 4
%                    symbols_tmp = symbols_tmp*2;  % rank 4: two layers for codeword 1 and 2
%                end
               switch obj.UEs(user_ind(uu)).traffic_model.type
                   case 'fullbuffer'
                       utility(uu,:,:) = c(uu,:,:).*symbols_tmp*max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp;
                   case 'voip'
                       utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,obj.UEs(user_ind(uu)).traffic_model.get_buffer_length)*-obj.lambda(user_ind(uu));
                   case 'ftp'
                       utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,obj.UEs(user_ind(uu)).traffic_model.get_buffer_length)*max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp;
                   case 'http'
                       utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,obj.UEs(user_ind(uu)).traffic_model.get_buffer_length)*max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp;
                   case 'video'
                       utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,obj.UEs(user_ind(uu)).traffic_model.get_buffer_length)*(max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp-obj.lambda(user_ind(uu)));
                   case 'gaming'
                       utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,obj.UEs(user_ind(uu)).traffic_model.get_buffer_length)*-obj.lambda(user_ind(uu));
                   otherwise
                       error('Traffic model not supported');
               end
           end
           RB_set = true(2,obj.RB_grid_size);
           RB_UEs = false(N_UE,2*obj.RB_grid_size,2);
           for rr = 1:2*obj.RB_grid_size % loop to assign resources             
               RB_set_temp = RB_set;
               tmp_c = zeros(2,obj.RB_grid_size);
               UE_ind = zeros(2,obj.RB_grid_size,2);
               for r_ = 1:sum(RB_set(:))
                   [r2,r1] = find(RB_set_temp,1,'first'); % currently considered resource
                   RB_set_temp(r2,r1) = false;
                       for rank_loop = min_rank:max_rank
                           UE_RI_ind = logical(RI_tmp == rank_loop);
                           if rank_loop == 1 % single stream transmission
                               [tmp_tmp_c,UE_ind_tmp] = max(utility(UE_RI_ind,r1+obj.RB_grid_size*(r2-1),1));
                               if tmp_tmp_c > tmp_c(r2,r1)
                                   tmp_c(r2,r1) = tmp_tmp_c;
                                   if tmp_c(r2,r1) ~= 10^-10 % no user is scheduled if tmp_c = 10^-10
                                        UE_ind_tmp = find(cumsum(UE_RI_ind) == UE_ind_tmp,1,'first');
                                        UE_ind(r2,r1,:) = [UE_ind_tmp,0];
                                   else
                                        UE_ind(r2,r1,:) = [0,0];
                                   end
                               end
                           elseif obj.MUMIMO % MUMIMO
                               pmi_tmp = squeeze(PMI_tmp(r1,r2,UE_RI_ind));
                               c_tmp = reshape(utility(UE_RI_ind,r1+obj.RB_grid_size*(r2-1),1:min(rank_loop,2)),sum(UE_RI_ind),[]);
                               for p_i = min(pmi_tmp):max(pmi_tmp)
                                   UE_PMI_ind = pmi_tmp == p_i;
                                   [c1,cw1] = max(c_tmp(UE_PMI_ind,1));
                                   [c2,cw2] = max(c_tmp(UE_PMI_ind,2));
                                   [calt,cwalt] = min(abs(c_tmp(UE_PMI_ind,2)-c2)); % alternative to c2 (if the buffer gets empty after assignment of the first stream)
                                   if c1+c2 > tmp_c(r2,r1)
                                       if c1 ~= 10^-10
                                            UE_ind_tmp1 = find(cumsum(UE_PMI_ind) == cw1,1,'first'); 
                                            UE_ind_tmp1 = find(cumsum(UE_RI_ind) == UE_ind_tmp1,1,'first');
                                       else
                                           UE_ind_tmp1 = 0;
                                       end
                                       if c2 ~= 10^-10
                                            UE_ind_tmp2 = find(cumsum(UE_PMI_ind) == cw2,1,'first');                                   
                                            UE_ind_tmp2 = find(cumsum(UE_RI_ind) == UE_ind_tmp2,1,'first');
                                       else
                                           UE_ind_tmp2 = 0;
                                       end
                                       if UE_ind_tmp1 ~= 0
                                       if UE_ind_tmp1 == UE_ind_tmp2 && ~strcmp(obj.UEs(user_ind(UE_ind_tmp1)).traffic_model.type,'fullbuffer') % if both streams go to a single user we have to check whether this really makes sense or the user buffer is too small
                                           bitnr = c(UE_ind_tmp1,r_,1)*(obj.Ns_RB-overhead(r_));
                                           symbols_tmp = symbols;
                                           if RI_tmp(UE_ind_tmp1) == 3
                                               symbols_tmp(:,:,1) = symbols_tmp(:,:,1)*2; % rank 3: two layers for codeword 1
                                               bitnr = bitnr*2;
                                           elseif RI_tmp(UE_ind_tmp1) == 4
                                               symbols_tmp = symbols_tmp*2;  % rank 4: two layers for codeword 1 and 2
                                               bitnr = bitnr*2;
                                           end
                                           switch obj.UEs(user_ind(UE_ind_tmp1)).traffic_model.type
%                                                case 'fullbuffer'
%                                                    util_tmp = c(uu,:,:).*symbols*max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp;
                                               case 'voip'
                                                   util_tmp = min(c(UE_ind_tmp1,r_,2).*symbols_tmp(1,r_,2),obj.UEs(user_ind(UE_ind_tmp1)).traffic_model.get_buffer_length-bitnr)*-obj.lambda(user_ind(UE_ind_tmp1));
                                               case 'ftp'
                                                   util_tmp = min(c(UE_ind_tmp1,r_,2).*symbols_tmp(1,r_,2),obj.UEs(user_ind(UE_ind_tmp1)).traffic_model.get_buffer_length-bitnr)*max(obj.av_throughput(user_ind(UE_ind_tmp1)),eps)^-alpha_tmp;
                                               case 'http'
                                                   util_tmp = min(c(UE_ind_tmp1,r_,2).*symbols_tmp(1,r_,2),obj.UEs(user_ind(UE_ind_tmp1)).traffic_model.get_buffer_length-bitnr)*max(obj.av_throughput(user_ind(UE_ind_tmp1)),eps)^-alpha_tmp;
                                               case 'video'
                                                   util_tmp = min(c(UE_ind_tmp1,r_,2).*symbols_tmp(1,r_,2),obj.UEs(user_ind(UE_ind_tmp1)).traffic_model.get_buffer_length-bitnr)*(max(obj.av_throughput(user_ind(UE_ind_tmp1)),eps)^-alpha_tmp-obj.lambda(user_ind(UE_ind_tmp1)));
                                               case 'gaming'
                                                   util_tmp = min(c(UE_ind_tmp1,r_,2).*symbols_tmp(1,r_,2),obj.UEs(user_ind(UE_ind_tmp1)).traffic_model.get_buffer_length-bitnr)*-obj.lambda(user_ind(UE_ind_tmp1));
                                               otherwise
                                                   error('Traffic model not supported');
                                           end
                                           if util_tmp < calt
                                               c2 = calt;
                                               cw2 = cwalt;
                                               if c2 ~= 10^-10
                                                    UE_ind_tmp2 = find(cumsum(UE_PMI_ind) == cw2,1,'first');                                   
                                                    UE_ind_tmp2 = find(cumsum(UE_RI_ind) == UE_ind_tmp2,1,'first');
                                               else
                                                    UE_ind_tmp2 = 0;
                                               end
                                           end
                                       end
                                       end
                                       tmp_c(r2,r1) = c1+c2;
                                       UE_ind(r2,r1,:) = [UE_ind_tmp1,UE_ind_tmp2];
                                   end
                                end
                           else % SUMIMO
                               [tmp_tmp_c,UE_ind_tmp] = max(sum(utility(UE_RI_ind,r1+obj.RB_grid_size*(r2-1),:),3));
                               if tmp_tmp_c > tmp_c(r2,r1)
                                   tmp_c(r2,r1) = tmp_tmp_c;
                                   if tmp_c(r2,r1) ~= 10^-10 % no user served
                                        UE_ind_tmp = find(cumsum(UE_RI_ind) == UE_ind_tmp,1,'first');
                                        UE_ind(r2,r1,:) = [UE_ind_tmp,UE_ind_tmp];
                                   else
                                        UE_ind(r2,r1,:) = [0,0];
                                   end
                               end 
                           end
                       end          
               end
               tmp_c = tmp_c';
               [~,max_ind] = max(tmp_c(:));
               [max_ind1,max_ind2] = ind2sub(size(tmp_c),max_ind);
               RB_set(max_ind2,max_ind1) = false;
               if UE_ind(max_ind2,max_ind1,1) ~= 0 % assign the resource for codeword 1 to the specific user
                   temp_ue = UE_ind(max_ind2,max_ind1,1);
                   mult = 1;
                   if RI_tmp(temp_ue) == 3 || RI_tmp(temp_ue) == 4
                       mult = 2;
                   end
                   RB_UEs(temp_ue,max_ind,1) = true;
                   if ~strcmp(obj.UEs(user_ind(temp_ue)).traffic_model.type,'fullbuffer')
                       if isinf(obj.UEs(user_ind(temp_ue)).traffic_model.bit_count)
                           obj.UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_ue,max_ind,1)*(obj.Ns_RB-overhead(max_ind))*mult,1);  % estimate the number of bits still left in the buffer
                       else
                           obj.UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_ue,max_ind,1)*(obj.Ns_RB-overhead(max_ind))*mult);
                       end
                   end
               end
               if UE_ind(max_ind2,max_ind1,2) ~= 0 % assign the resource for codeword 2 to the specific user
                   temp_ue = UE_ind(max_ind2,max_ind1,2);
                   mult = 1;
                   if RI_tmp(temp_ue) == 4
                       mult = 2;
                   end
                   RB_UEs(temp_ue,max_ind,2) = true;
                   if ~strcmp(obj.UEs(user_ind(temp_ue)).traffic_model.type,'fullbuffer')
                       if isinf(obj.UEs(user_ind(temp_ue)).traffic_model.bit_count)
                           obj.UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_ue,max_ind,2)*(obj.Ns_RB-overhead(max_ind))*mult,1);  % estimate the number of bits still left in the buffer
                       else
                           obj.UEs(user_ind(temp_ue)).traffic_model.coarse_decrease(c(temp_ue,max_ind,2)*(obj.Ns_RB-overhead(max_ind))*mult);
                       end
                   end
%                    RB_UEs(UE_ind(max_ind2,max_ind1,2),max_ind,2) = true;
               end
               for uu = 1:N_UE % update user utility for all resources (necessary if user buffer gets empty during resource allocation
                   if ~strcmp(obj.UEs(user_ind(uu)).traffic_model.type,'fullbuffer')
                       bitnr = obj.UEs(user_ind(uu)).traffic_model.bit_count; % number of bits left in the buffer
                       if  bitnr <= 0 % this user does not have any data left
                           utility(uu,:,:) = 10^-10;
                       else
                           symbols_tmp = symbols;
                           if RI_tmp(uu) == 3
                               symbols_tmp(:,:,1) = symbols_tmp(:,:,1)*2; % rank 3: two layers for codeword 1
                           elseif RI_tmp(uu) == 4
                               symbols_tmp = symbols_tmp*2;  % rank 4: two layers for codeword 1 and 2
                           end
                           switch obj.UEs(user_ind(uu)).traffic_model.type
                               case 'fullbuffer' % there is nothing to do for fullbuffer simulations (the utilities do not change)
                               case 'voip'
                                   utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,bitnr)*-obj.lambda(user_ind(uu));
                               case 'ftp'
                                   utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,bitnr)*max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp;
                               case 'http'
                                   utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,bitnr)*max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp;
                               case 'video'
                                   utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,bitnr)*(max(obj.av_throughput(user_ind(uu)),eps)^-alpha_tmp-obj.lambda(user_ind(uu)));
                               case 'gaming'
                                   utility(uu,:,:) = min(c(uu,:,:).*symbols_tmp,bitnr)*-obj.lambda(user_ind(uu));
                               otherwise
                                   error('Traffic model not supported');
                           end
                       end
                   end
               end
           end
           RBs(:,1) = reshape(RB_UEs(:,:,1),[],1);
           RBs(:,2) = reshape(RB_UEs(:,:,2),[],1);
           if ~sum(RBs(:,2)) % just one layer is used
               RBs = RBs(:,1);
           end
%            RBs = RBs(:,1);
       end
   end
end 

