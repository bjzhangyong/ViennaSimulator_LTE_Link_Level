classdef fixedScheduler < network_elements.lteScheduler
% Scheduler that assigns a fixed and predefined amount of RBs to each user.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       CQI_mapping_data
       SINR_averager
   end

   methods
       function obj = fixedScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data)
           
           % Fill in basic parameters (handled by the superclass constructor)
           obj = obj@network_elements.lteScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params);

           if length(scheduler_params.UE_assignment)~=length(UEs_to_be_scheduled)
               error('You must define the number of assigned RBs for EACH user. %d UEs were assigned, but %d exist',...
                   length(scheduler_params.UE_assignment),length(UEs_to_be_scheduled));
           elseif sum(scheduler_params.UE_assignment)>RB_grid_size*2
               error('You cannot assign more RBs that the ones available. Assigned %d RBs. Only %d available',...
                   sum(scheduler_params.UE_assignment),RB_grid_size*2);
           end

           switch scheduler_params.assignment
               case 'static'
                   % Static CQI assignment
                   obj.static_scheduler = true;
               otherwise
                   obj.static_scheduler = false;
           end

           % Get a vector of scheduling params (one for each UE)
           % initialized to the values that we want
           obj.UE_static_params = obj.get_initialized_UE_params(scheduler_params,CQI_params);
           obj.SINR_averager = averager;
           obj.CQI_mapping_data = mapping_data;
           % Fill in the RB allocation grid for each user (and codeword)
           %UE_mapping_all_UEs = zeros(RB_grid_size,2,obj.maxCodewords);
           UE_mapping_all_UEs = zeros(RB_grid_size,2);
           starting_rb = 1;
           for u_=1:obj.nUEs
               for rb_=starting_rb:starting_rb+scheduler_params.UE_assignment(u_)-1
                   % NOTE: Same RB assignment for both codewords.
                   UE_mapping_all_UEs(rb_) = u_;
%                    for cw_=2:scheduler_params.nCodewords(u_)
%                        UE_mapping_all_UEs(:,cw_) = UE_mapping_all_UEs(:,1);
%                    end
               end
               starting_rb = starting_rb+scheduler_params.UE_assignment(u_);
           end
           % Assign the static scheduling parameters for each user
           for u_=1:obj.nUEs
               obj.UE_static_params(u_).UE_mapping = (UE_mapping_all_UEs==u_);
               obj.UE_static_params(u_).assigned_RBs = squeeze(sum(sum(obj.UE_static_params(u_).UE_mapping,1),2));
           end
       end
       
       % Generate UE scheduling.
       %
       % NOTE: if you modify the returned
       % UE_scheduling objects and the scheduler happens to be static, you
       % will modify the static scheduling info stored in the scheduler
       % permanently!! To avoid this I should have to copy the scheduler
       % info every time, which I don't want to, as I don't see any reason
       % to change key parameters of the static scheduler during a static
       % scheduling simulation. Another option is to have the scheduling
       % info by value instead of by reference (handle), which is even a
       % worse solution.
       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           
           UE_scheduling = obj.UE_static_params;
           for uu = 1:size(UE_output,2)
               if ~obj.static_scheduler
                   UE_scheduling(uu).CQI_params = [];
               end
               if ~isempty(UE_output(uu).PMI)   % check wheter PMI is fed back
                   if ~isempty(UE_output(uu).RI)    % check wheter RI is fed back
                        UE_scheduling(uu).nLayers = UE_output(uu).RI;
                        UE_scheduling(uu).nCodewords = min(2,UE_output(uu).RI);
                   end
               UE_scheduling(uu).PMI = UE_output(uu).PMI;
               end
%                if ~isempty(UE_output(uu).CQI)
                   UE_scheduling(uu).cqi = zeros(1,UE_scheduling(uu).nCodewords);
                   for i1 = 1:UE_scheduling(uu).nCodewords
                       if ~obj.static_scheduler
                           CQI_all_temp = UE_output(uu).CQI(:,:,i1);
                           [val,~] = max(abs(diff(CQI_all_temp(:))));
                           if val ~= 0 && sum(UE_scheduling(uu).UE_mapping(:)) ~= 0 % different CQI values at different RBs --> averaging necessary!
                               CQI_temp = CQI_all_temp(UE_scheduling(uu).UE_mapping);
                               SINRs = zeros(size(obj.CQI_mapping_data));
                               %                             SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data(mod(CQI_temp,20)+1))/10),0:15);
                               SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data(mod(CQI_temp,20)+1)+obj.CQI_mapping_data(mod(CQI_temp,20)+2))/20),0:15); % this version is less conservative
                               SINRs = SINR_temp(:);
                               temp = zeros(size(SINRs));
                               temp(obj.CQI_mapping_data(1:16) <= SINRs) = 1;
                               temp_CQI = find(temp,1,'last')-1;
                               if temp_CQI
                                   UE_scheduling(uu).cqi(i1) = temp_CQI;
                               else
                                   UE_scheduling(uu).cqi(i1) = 20; % this is the rate 0 CQI
                               end
                           else
                               UE_scheduling(uu).cqi(i1) = squeeze(UE_output(uu).CQI(1,1,i1));
                           end
                           if logical(~UE_scheduling(uu).cqi(i1)) || (UE_scheduling(uu).cqi(i1) == 20)
                               UE_scheduling(uu).cqi(i1) = 1;
                           end
                       else
                           UE_scheduling(uu).cqi(i1) = UE_scheduling(uu).CQI_params.CQI;
                       end
                       UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi(i1),obj.CQI_params)];
                   end
%                else
%                    for i1 = 1:UE_scheduling(uu).nCodewords
%                         UE_scheduling(uu).cqi(i1) = obj.cqi_i;
%                         UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi(i1),obj.CQI_params)];
%                    end
%                end
           end
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
       end
   end
end 
