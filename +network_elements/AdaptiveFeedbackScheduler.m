classdef AdaptiveFeedbackScheduler < network_elements.lteScheduler
% Scheduler that assigns a fixed and predefined amount of RBs to each user
% and dynamically adjusts PMI and RI
% Stefan Schwarz 
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       SINR_averager
       CQI_mapping_data
       RB_gridsize
       cqi_i
       alphabets
   end

   methods
       function obj = AdaptiveFeedbackScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params,averager,mapping_data,varargin)
           % Fill in basic parameters (handled by the superclass constructor)
           obj = obj@network_elements.lteScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params);
           obj.static_scheduler = true;
           obj.SINR_averager = averager;
           obj.CQI_mapping_data = mapping_data;
           % Get a vector of scheduling params (one for each UE)
           % initialized to the values that we want
           obj.UE_static_params = obj.get_initialized_UE_params(scheduler_params,CQI_params);
           obj.cqi_i = scheduler_params.cqi;
           obj.RB_gridsize = RB_grid_size;
           if ~isempty(varargin)
               obj.alphabets = varargin{1};
           else
               obj.alphabets = [];
           end
       end
       
       function RB_allocation(obj)
           % Fill in the RB allocation grid for each user (and codeword)
           RB_grid_size = obj.RB_gridsize;
           UE_mapping_all_UEs = zeros(2,RB_grid_size);  % the two corresponds to the number of slots; keep it that way to be fair, otherwise UEs do not see the same amount of sync signals
           RB_num_high = ceil(RB_grid_size*2/obj.nUEs);
           RB_num_low = floor(RB_grid_size*2/obj.nUEs);
           if RB_num_high ~= RB_num_low
               change = (obj.nUEs*RB_num_high-RB_grid_size*2)/(RB_num_high-RB_num_low);
           else 
               change = obj.nUEs+1;
           end
           rand_users = randperm(obj.nUEs); % randomize user allocation
%            rand_users = 1:obj.nUEs;
           for i1 = 1:obj.nUEs
              if i1 <= change
                  UE_mapping_all_UEs((i1-1)*RB_num_low+1:i1*RB_num_low) = rand_users(i1);
              else
                  UE_mapping_all_UEs(change*RB_num_low+(i1-change-1)*RB_num_high+1:change*RB_num_low+(i1-change)*RB_num_high) = rand_users(i1);
              end
           end
           UE_mapping_all_UEs = UE_mapping_all_UEs.';
           % Assign the static scheduling parameters for each user
           for u_=1:obj.nUEs
               obj.UE_static_params(u_).UE_mapping = (UE_mapping_all_UEs==u_);
               obj.UE_static_params(u_).assigned_RBs = squeeze(sum(sum(obj.UE_static_params(u_).UE_mapping,1),2));
           end
       end

       function UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data,cell_genie,PBCHsyms)
           % this function schedules the user on his statically configured
           % resource blocks and dynamically adjusts the transmission rank and
           % precoding matrix indicator
%            UE_scheduling = obj.UE_static_params;

           obj.RB_allocation; % redo the resource allocation everytime to average the number of resources a user gets (if N_RB/N_UE is noninteger)

           UE_scheduling = obj.UE_static_params;
           N_UE = size(UE_output,2);
           N_RB = size(UE_output(1).CQI,1)*2;
           UE_scheduling = obj.UE_static_params;
           RBs = zeros(N_UE*N_RB,1);
           for uu = 1:N_UE
               UE_scheduling(uu).CQI_params = [];
               if ~isempty(UE_output(uu).PMI)   % check wheter PMI is fed back
                   if ~isempty(UE_output(uu).RI)    % check wheter RI is fed back
                        UE_scheduling(uu).nLayers = UE_output(uu).RI;
                        UE_scheduling(uu).nCodewords = min(2,UE_output(uu).RI);
                   end
               UE_scheduling(uu).PMI = UE_output(uu).PMI;
               end
               RBs(uu:N_UE:end) = UE_scheduling(uu).UE_mapping(:);
           end
           UE_scheduling = obj.set_cqi(UE_scheduling,1:N_UE,RBs(1:end),N_UE,N_RB,UE_output);
%            for uu = 1:size(UE_output,2)
%                UE_scheduling(uu).CQI_params = [];
%                if ~isempty(UE_output(uu).PMI)   % check wheter PMI is fed back
%                     UE_scheduling(uu).PMI = UE_output(uu).PMI;
%                end
%                if ~isempty(UE_output(uu).RI)    % check wheter RI is fed back
%                     UE_scheduling(uu).nLayers = UE_output(uu).RI;
%                     UE_scheduling(uu).nCodewords = min(2,UE_output(uu).RI);
%                end
% %                if ~isempty(UE_output(uu).CQI)
%                    UE_scheduling(uu).cqi = zeros(1,UE_scheduling(uu).nCodewords);
%                    for i1 = 1:UE_scheduling(uu).nCodewords
%                        CQI_all_temp = UE_output(uu).CQI(:,:,i1);
%                        [val,~] = max(abs(diff(CQI_all_temp(:))));
%                        if val ~= 0 && sum(UE_scheduling(uu).UE_mapping(:)) ~= 0 % different CQI values at different RBs --> averaging necessary!
%                             CQI_temp = CQI_all_temp(UE_scheduling(uu).UE_mapping);
%                             SINRs = zeros(size(obj.CQI_mapping_data));
% %                             SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data(mod(CQI_temp,20)+1))/10),0:15);
%                             SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data(mod(CQI_temp,20)+1)+obj.CQI_mapping_data(mod(CQI_temp,20)+2))/20),0:15); % this version is less conservative
%                             SINRs = SINR_temp(:);
%                             temp = zeros(size(SINRs));
%                             temp(obj.CQI_mapping_data(1:16) <= SINRs) = 1;
%                             temp_CQI = find(temp,1,'last')-1;
%                             if temp_CQI
%                                 UE_scheduling(uu).cqi(i1) = temp_CQI;
%                             else
%                                 UE_scheduling(uu).cqi(i1) = 20; % this is the rate 0 CQI
%                             end
%                        else
%                            UE_scheduling(uu).cqi(i1) = squeeze(UE_output(uu).CQI(1,1,i1));
%                        end
%                        if logical(~UE_scheduling(uu).cqi(i1)) || (UE_scheduling(uu).cqi(i1) == 20)
%                             UE_scheduling(uu).cqi(i1) = 1;
%                        end
%                        UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi(i1),obj.CQI_params)];
%                    end
% %                else
% %                    for i1 = 1:UE_scheduling(uu).nCodewords
% %                         UE_scheduling(uu).cqi(i1) = obj.cqi_i;
% %                         UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi(i1),obj.CQI_params)];
% %                    end
% %                end
%            end
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
       end
   end
end 
