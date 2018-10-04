classdef lteScheduler < handle
% Implements methods common to all the schedulers.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       maxCodewords = 2;  % Maximum number of codewords that can be transmitted per TTI. For LTE, that's 2
       static_scheduler   % Whether this scheduler is static or dynamic
       UEs                % List of UEs to schedule
       nUEs               % Total number of UEs -> length(UEs)
       UE_static_params   % In case this is a static scheduler, this variable stores for each
                          % UE the scheduling parameters that will be user.
                          % The only parameter that cannot be statically
       Ns_RB              % Number of symbols in one RB (needed to calculate the TB sizes)
       max_HARQ_retx      % Max num of HARQ retransmissions, NOT including the original tx. 0, 1, 2 or 3
       attached_eNodeB    % eNodeB to which this scheduler is attached
       RB_grid_size       % Size of the Resource Block grid
       
       CQI_params         % CQI parameters for all possible MCSs
       
       zero_delay         % Specify whether there is zero delay for the CQI reports. Note that zero-delay
                          % CQI reports are always calculated with perfect
                          % channel knowledge.
                          
       CQI_mapping_params % parameters needed to perform the CQI mapping
       
       UE_specific        % direct reference to the HARQ processes from the eNodeB
       first              % indicator variable
       av_const           % averaging constant for the exponential throughput averaging
       overhead_ref       % overhead due to reference symbols per RB
       overhead_sync      % overhead due to synchronization symbols per RB
       RB_count
   end

   methods
       % Superclass constructor which is called by all subclasses
       function obj = lteScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params)
           
           obj.RB_grid_size  = RB_grid_size;
           obj.Ns_RB         = Ns_RB;
           obj.UEs           = UEs_to_be_scheduled;
           obj.nUEs          = length(UEs_to_be_scheduled);
           
           obj.max_HARQ_retx = scheduler_params.max_HARQ_retx;
           obj.zero_delay    = scheduler_params.zero_delay;
           obj.CQI_params    = CQI_params;
           
           obj.CQI_mapping_params = scheduler_params.CQI_mapping_params;
           obj.UE_specific        = scheduler_params.UE_specific;
           obj.av_const = scheduler_params.av_window;
           obj.overhead_ref = scheduler_params.overhead_ref;
           obj.overhead_sync = scheduler_params.overhead_sync;
           obj.first = true;
           obj.RB_count = zeros(obj.nUEs,1);
       end
       
       % Returns nUEs empty scheduler parameter objects
       function empty_UE_scheduling_params = get_new_UE_params(obj)
           empty_UE_scheduling_params = network_elements.UE_scheduling_params;
           for u_ = 2:obj.nUEs
               empty_UE_scheduling_params(u_) = network_elements.UE_scheduling_params;
           end
       end
       
       % Returns nUEs scheduler parameter objects initialized to the given
       % scheduler params. Only the UE allocation is not initialized. This
       % function is useful to program new static schedulers (less code
       % needed there). If some parameters are not assigned in the
       % scheduler_params struct, they will not be initialized
       function initialized_UE_scheduling_params = get_initialized_UE_params(obj,scheduler_params,CQI_params)
           
           % Generate all of the static assignments as specified in scheduler_params
           
           % TX mode, number of layers and number of codewords
           if isfield(scheduler_params,'tx_mode')
               tx_mode = scheduler_params.tx_mode; % tx mode for each UE (vector)
               nLayers    = scheduler_params.nLayers;
               nCodewords = scheduler_params.nCodewords;
           else
               tx_mode = [];
               nLayers = [];
               nCodewords = [];
           end
           % CQI assignment (if set, assume the same for all users)
           if isfield(scheduler_params,'cqi') && strcmp(scheduler_params.assignment,'static')
               if length(scheduler_params.cqi) == 1
                   cqi = scheduler_params.cqi * ones(nCodewords(1),1); % CQI for each user (n codewords)
               elseif length(scheduler_params.cqi) == nCodewords(1)
                   cqi = scheduler_params.cqi;
               else
                   error('something is wrong with number of CQIs');
               end
               cqi_params = LTE_common_get_CQI_params(cqi,CQI_params);
%            elseif strcmp(scheduler_params.assignment,'semi static') || strcmp(scheduler_params.type,'proportional fair') || strcmp(scheduler_params.type,'best cqi')
%                cqi = 1 * ones(nCodewords(1),1); % CQI for each user (n codewords)
%                cqi_params = LTE_common_get_CQI_params(cqi,CQI_params);
           else
               cqi = [];
               cqi_params = [];
           end

           % PMI assignment
           if isfield(scheduler_params,'PMI')
               PMI = scheduler_params.PMI;
           else
               PMI = [];
           end
           
           if isfield(scheduler_params,'CDD')
               CDD = scheduler_params.CDD;
           else
               CDD = [];
           end
           
           initialized_UE_scheduling_params = obj.get_new_UE_params;
           
           % Assign the static scheduling parameters for each user
           for u_=1:obj.nUEs
               
               if ~isempty(tx_mode)
                   initialized_UE_scheduling_params(u_).tx_mode = tx_mode(u_);
                   initialized_UE_scheduling_params(u_).nLayers = nLayers(u_);
                   initialized_UE_scheduling_params(u_).nCodewords = nCodewords(u_);
               end
               
               if ~isempty(cqi)
                   initialized_UE_scheduling_params(u_).cqi = cqi;
                   initialized_UE_scheduling_params(u_).CQI_params = cqi_params;
               end
               
               initialized_UE_scheduling_params(u_).UE_mapping = [];
               
               % MIMO parameters
               if ~isempty(PMI)
                   initialized_UE_scheduling_params(u_).PMI = PMI;
               end
               if ~isempty(CDD)
                   initialized_UE_scheduling_params(u_).CDD = CDD;
               end
               % NOTE: Please note that due to ease of implementation
               % these parameters have been included here, but take
               % into account that these parameters are in NO WAY
               % static!!!
               initialized_UE_scheduling_params(u_).codebook_index = [];
               initialized_UE_scheduling_params(u_).PRE = [];
%                initialized_UE_scheduling_params(u_).W = [];
%                initialized_UE_scheduling_params(u_).U = [];
               initialized_UE_scheduling_params(u_).freq_indices = [];
           end
       end

       % This fills in the number of bits allocated to this user according
       % to the specified scheduling parameters and also updates the HARQ
       % process objects with the scheduling information
%        function calculate_allocated_bits(obj,UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements)
%            
%            % You can access the properties of the eNodeB like this
%            BS_nTX = obj.attached_eNodeB.nTX;
%            
%            for u_=1:length(UE_scheduling)
%                
%                % Calculation of the TB size
%                if(subframe_corr == 1 || subframe_corr == 6)
%                    %lenghts of primary and secondary synchronization channels (symbols)
%                    sync_symbols = sum(sum(SyncUsedElements(UE_scheduling(u_).UE_mapping)));
%                else
%                    sync_symbols = 0;
%                end
%               
%                if(UE_scheduling(u_).assigned_RBs)
% 
%                    switch UE_scheduling(u_).nCodewords
%                        case 1
%                            assigned_RBs = UE_scheduling(u_).assigned_RBs*2;
%                        case 2
%                            assigned_RBs = UE_scheduling(u_).assigned_RBs;
%                    end
% 
%                    UE_scheduling(u_).N_coded_bits = [UE_scheduling(u_).CQI_params.modulation_order] * (assigned_RBs * (obj.Ns_RB - total_no_refsym) - sync_symbols);
%                    UE_scheduling(u_).N_data_bits  = 8*round(1/8 * UE_scheduling(u_).N_coded_bits .* [UE_scheduling(u_).CQI_params.coding_rate_x_1024] / 1024)-24; % calculate G based on TB_size and the target rate
%                    
%                    if(UE_scheduling(u_).N_data_bits < 16)
%                        UE_scheduling(u_).N_coded_bits = 0;
%                        UE_scheduling(u_).N_data_bits  = 0;
%                        UE_scheduling(u_).UE_mapping   = false(obj.RB_grid_size,2);
%                        UE_scheduling(u_).assigned_RBs = 0;
%                        % clear parameters for the case when the user is not scheduled (when the number of data bits is less than 16)
%                        UE_scheduling(u_).cqi = 0 * ones(obj.maxCodewords,1);
%                        UE_scheduling(u_).CQI_params = [];
%                    end
%                else
%                    UE_scheduling(u_).N_coded_bits = 0;
%                    UE_scheduling(u_).N_data_bits = 0;
%                end
%  
%                % Adjust number of bits for the SM case with 4 antennas
%                if (UE_scheduling(u_).tx_mode ~= 1 && UE_scheduling(u_).tx_mode ~= 2 && BS_nTX == 4)
%                    switch UE_scheduling(u_).nLayers
%                        case 2
%                            if (UE_scheduling(u_).nCodewords == 1)
%                                UE_scheduling(u_).N_coded_bits = UE_scheduling(u_).N_coded_bits*2;
%                                UE_scheduling(u_).N_data_bits = UE_scheduling(u_).N_data_bits*2;
%                            end
%                        case 3
%                            UE_scheduling(u_).N_coded_bits(2) = UE_scheduling(u_).N_coded_bits(2)*2;
%                            UE_scheduling(u_).N_data_bits(2) = UE_scheduling(u_).N_data_bits(2)*2;
%                        case 4
%                            UE_scheduling(u_).N_coded_bits = UE_scheduling(u_).N_coded_bits*2;
%                            UE_scheduling(u_).N_data_bits = UE_scheduling(u_).N_data_bits*2;
%                    end
%                end
% 
%                % Adjust information from the HARQ processes
%                for cw_=1:UE_scheduling(u_).nCodewords
%                    % Update HARQ process information
%                    obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs = UE_scheduling(u_).assigned_RBs;
%                    obj.UE_specific(u_).current_HARQ_process(cw_).cqi          = UE_scheduling(u_).cqi;
%                    % rv_idx information of the HARQ process is updated outside of the scheduler
%                    % Update scheduler information only
%                    UE_scheduling(u_).HARQ_process_id(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).id;
%                    UE_scheduling(u_).rv_idx(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx;
%                end
%                
%                for cw_=UE_scheduling(u_).nCodewords+1:size(obj.UE_specific(u_).HARQ_processes,1)
%                    obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs = 0;
%                    obj.UE_specific(u_).current_HARQ_process(cw_).cqi          = 0;
%                    % rv_idx information of the HARQ process is updated outside of the scheduler
%                    % Update scheduler information only
%                    UE_scheduling(u_).rv_idx(cw_) = 0;
%                end
%                
%            end
%        end
       
       function calculate_allocated_bits(obj,UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,ChUsedElements)
           
           % You can access the properties of the eNodeB like this
           BS_nTX = obj.attached_eNodeB.nTX;
           for u_=1:length(UE_scheduling)
               CHsyms = sum(sum(ChUsedElements(UE_scheduling(u_).UE_mapping)));
               % Calculation of the TB size
               if(subframe_corr == 1 || subframe_corr == 6)
                   %lenghts of primary and secondary synchronization channels (symbols)
                   sync_symbols = sum(sum(SyncUsedElements(UE_scheduling(u_).UE_mapping)));               
               else
                   sync_symbols = 0;
               end

               if(UE_scheduling(u_).assigned_RBs)
                   assigned_RBs = UE_scheduling(u_).assigned_RBs;
%                    switch UE_scheduling(u_).nCodewords
%                        case 1
%                            assigned_RBs = UE_scheduling(u_).assigned_RBs*2;
%                        case 2
%                            assigned_RBs = UE_scheduling(u_).assigned_RBs;
%                    end
                   obj.RB_count(u_) = obj.RB_count(u_)+assigned_RBs;
                   UE_scheduling(u_).N_coded_bits = [UE_scheduling(u_).CQI_params.modulation_order] * (assigned_RBs * (obj.Ns_RB - total_no_refsym) - sync_symbols - CHsyms);
                   UE_scheduling(u_).N_data_bits  = 8*round(1/8 * UE_scheduling(u_).N_coded_bits .* [UE_scheduling(u_).CQI_params.coding_rate_x_1024] / 1024)-24; % calculate G based on TB_size and the target rate
                   if assigned_RBs * (obj.Ns_RB - total_no_refsym) - sync_symbols - CHsyms == 0
                       UE_scheduling(u_).assigned_RBs = zeros(size(assigned_RBs));
                       UE_scheduling(u_).N_data_bits = 0;
                       UE_scheduling(u_).N_coded_bits = 0;
                   else
                       for i_ = 1:length(UE_scheduling(u_).N_data_bits)
                           while 1
                                if(UE_scheduling(u_).N_data_bits(i_) < 16)
                                    UE_scheduling(u_).cqi(i_) = UE_scheduling(u_).cqi(i_)+1;
                                    UE_scheduling(u_).CQI_params(i_) = LTE_common_get_CQI_params(UE_scheduling(u_).cqi(i_),obj.CQI_params);
                                    UE_scheduling(u_).N_coded_bits(i_) = UE_scheduling(u_).CQI_params(i_).modulation_order * (assigned_RBs * (obj.Ns_RB - total_no_refsym) - sync_symbols - CHsyms);
                                    UE_scheduling(u_).N_data_bits(i_)  = 8*round(1/8 * UE_scheduling(u_).N_coded_bits(i_) .* UE_scheduling(u_).CQI_params(i_).coding_rate_x_1024 / 1024)-24; % calculate G based on TB_size and the target rate
                                else
                                    break
                                end
                           end
                       end
                   end
                   if sum(UE_scheduling(u_).N_data_bits) == 0
                       UE_scheduling(u_).assigned_RBs = 0;
                       UE_scheduling(u_).UE_mapping   = false(obj.RB_grid_size,2);
                   end
               else
                   UE_scheduling(u_).N_coded_bits = 0;
                   UE_scheduling(u_).N_data_bits = 0;
               end

               % Adjust number of bits for the SM case with 4 antennas
               if (UE_scheduling(u_).tx_mode ~= 1 && UE_scheduling(u_).tx_mode ~= 2 && BS_nTX == 4)
                   switch UE_scheduling(u_).nLayers
                       case 2
                           if (UE_scheduling(u_).nCodewords == 1)
                               UE_scheduling(u_).N_coded_bits = UE_scheduling(u_).N_coded_bits*2;
                               UE_scheduling(u_).N_data_bits = UE_scheduling(u_).N_data_bits*2;
                           end
                       case 3
                           UE_scheduling(u_).N_coded_bits(2) = UE_scheduling(u_).N_coded_bits(2)*2;
                           UE_scheduling(u_).N_data_bits(2) = UE_scheduling(u_).N_data_bits(2)*2;
                       case 4
                           UE_scheduling(u_).N_coded_bits = UE_scheduling(u_).N_coded_bits*2;
                           UE_scheduling(u_).N_data_bits = UE_scheduling(u_).N_data_bits*2;
                   end
               end

               % Adjust information from the HARQ processes
               UE_scheduling(u_).N_used_bits = zeros(UE_scheduling(u_).nCodewords,1); 
               for cw_=1:UE_scheduling(u_).nCodewords
                   % Update HARQ process information
                   obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs = UE_scheduling(u_).assigned_RBs;
                   obj.UE_specific(u_).current_HARQ_process(cw_).cqi          = UE_scheduling(u_).cqi;
                   % rv_idx information of the HARQ process is updated outside of the scheduler
                   % Update scheduler information only
                   UE_scheduling(u_).HARQ_process_id(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).id;
                   UE_scheduling(u_).rv_idx(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx;
                   if UE_scheduling(u_).N_data_bits(cw_) ~= 0
                       if strcmp(obj.UEs(u_).traffic_model.type,'voip') || strcmp(obj.UEs(u_).traffic_model.type,'video') || strcmp(obj.UEs(u_).traffic_model.type,'gaming')
                            packet_parts = obj.UEs(u_).traffic_model.decrease_packets(UE_scheduling(u_).N_data_bits(cw_));
                            if ~isempty(packet_parts)
                                UE_scheduling(u_).N_used_bits(cw_) = sum(packet_parts.get_size);
                                obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts = packet_parts;
                            else
                                UE_scheduling(u_).N_used_bits(cw_) = 0;
                            end
                       else
                           obj.UEs(u_).traffic_model.decrease_packets(UE_scheduling(u_).N_data_bits(cw_));
                           UE_scheduling(u_).N_used_bits(cw_) = min(obj.UEs(u_).traffic_model.get_buffer_length,UE_scheduling(u_).N_data_bits(cw_));
                       end
                   else
                       UE_scheduling(u_).N_used_bits(cw_) =  0;
                   end
               end
               
               for cw_=UE_scheduling(u_).nCodewords+1:size(obj.UE_specific(u_).HARQ_processes,1)
                   obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs = 0;
                   obj.UE_specific(u_).current_HARQ_process(cw_).cqi          = 0;
                   % rv_idx information of the HARQ process is updated outside of the scheduler
                   % Update scheduler information only
                   UE_scheduling(u_).rv_idx(cw_) = 0;
               end
           end
       end

       function feedback_users_cqi_vec_all = get_all_UE_feedback(obj,UE_output)
           % Get feedback from all users. If we are using zero delay, then
           % calculate the CQIs from the genie SNRs. If not, use the
           % received UE feedback
           feedback_users_cqi_vec_all = zeros(obj.nUEs,obj.RB_grid_size);
           if obj.zero_delay
               % TODO: get genie information and calculate CQI
           else
               for u_=1:obj.nUEs
                   current_UE_feedback = UE_output(u_).CQI_feedback;
                   % It could be that there is no feedback (eg. Not yet
                   % arrived). Then set CQIs to all 0 (no transmission).
                   if isempty(current_UE_feedback)
                       feedback_users_cqi_vec_all(u_,:) = 0;
                   else
                       feedback_users_cqi_vec_all(u_,:) = current_UE_feedback;
                   end
               end
           end
       end
       
       function [TP,varargout] = compute_av_throughput(obj,UE_output,TP_old,uu,N_UE,varargin)
               TP = 0;
               if ~isempty(UE_output.rx_data_bits)
                   if ~obj.first
                       for cc = 1:length(UE_output.rx_data_bits)
                            TP = TP + UE_output.ACK(cc)*length(UE_output.rx_data_bits{cc});
                       end
                   else
                       for cc = 1:length(UE_output.rx_data_bits)
                            TP = UE_output.ACK(cc)*length(UE_output.rx_data_bits{cc});
                       end
                      if uu == N_UE
                            obj.first = false;
                      end
                   end   
               end
%               if sum(TP_old)
                  if ~isempty(varargin)
                      varargout{1} = [varargin{1}(2:end),TP];
    %                    varargout{1} = (varargin{1}*varargin{2}+ TP)/(varargin{2}+1);
        %                      varargout{1} = TP;
                  end
%               else
%                   varargout{1} = TP;
%               end
              if sum(TP_old)
                  TP = TP_old*(1-1/obj.av_const)+TP*1/obj.av_const; 
              end
       end
       
       function UE_scheduling = set_cqi(obj,UE_scheduling,user_ind,RBs,N_UE,N_RB,UE_output)
          for uu = 1:N_UE
              temp_ind = find(user_ind == uu);
              temp_mapping = full(RBs(temp_ind:N_UE:end));
              UE_scheduling(uu).UE_mapping = logical(reshape(temp_mapping,N_RB/2,2));
              UE_scheduling(uu).assigned_RBs = squeeze(sum(sum(obj.UE_static_params(uu).UE_mapping,1),2));
              UE_scheduling(uu).cqi = zeros(1,UE_scheduling(uu).nCodewords);
              if UE_scheduling(uu).assigned_RBs ~= 0
                  for i1 = 1:UE_scheduling(uu).nCodewords
                        CQI_all_temp = UE_output(uu).CQI(:,:,i1);
                        [val,~] = max(abs(diff(CQI_all_temp(:))));
                        if val % different CQI values at different RBs --> averaging necessary!
                            CQI_temp = CQI_all_temp(UE_scheduling(uu).UE_mapping);
                            SINRs = zeros(size(obj.CQI_mapping_data));
    %                                     SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data.table(mod(CQI_temp,20)+1))/10),0:15,obj.alphabets);
    %                         SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data.table(mod(CQI_temp,20)+1)+obj.CQI_mapping_data.table(mod(CQI_temp,20)+2))/20),0:15,obj.alphabets);
                            SINR_temp = obj.SINR_averager.average(10.^((obj.CQI_mapping_data.table(mod(CQI_temp,20)+1)+(obj.CQI_mapping_data.table(mod(CQI_temp,20)+2)-obj.CQI_mapping_data.table(mod(CQI_temp,20)+1))/2.5)/10),0:15,obj.alphabets); % this version is less conservative
                            SINRs = SINR_temp(:);
                            temp = zeros(size(SINRs));
                            temp(obj.CQI_mapping_data.table(1:16) <= SINRs) = 1;
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
                        UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi(i1),obj.CQI_params)];
                  end
              else
                  UE_scheduling(uu).cqi(1) = 1;
                  UE_scheduling(uu).CQI_params = [UE_scheduling(uu).CQI_params,LTE_common_get_CQI_params(UE_scheduling(uu).cqi(1),obj.CQI_params)];
              end
           end
       end
       
       function [UE_scheduling,c,user_ind,varargout] = set_pmi_ri(obj,UE_scheduling,N_UE,N_RB,UE_output)
           a = zeros(N_UE,1);
           b = zeros(N_UE,1);
           c = zeros(N_UE,N_RB);
           temp_UE = struct([]);
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
                a(uu) = obj.SINR_averager.MI_data(1).k(CQI_bar);
                b(uu) = obj.SINR_averager.MI_data(1).d(CQI_bar);
           end
           varargout{1} = a;
           varargout{2} = b;
           user_ind = randperm(N_UE);
%            user_ind = 1:N_UE;
           for rb = 1:N_RB
               c_count = 0;
               for uu = user_ind
                    c_count = c_count+1;
                    c(c_count,rb) = a(uu)* temp_UE(uu).CQI(rb)+b(uu); % efficiency in bit/channel use
               end
           end
%            c_count = 0;
%            for u_ = user_ind% calculate sum efficiencies for both codewords on every RB
%                c_count = c_count+1;
%                UE_scheduling(u_).UE_mapping = false(obj.RB_grid_size,2);
%                c(c_count,:) = sum(reshape([obj.CQI_params(UE_output(u_).CQI(:,:,1:UE_scheduling(u_).nCodewords)).efficiency],obj.RB_grid_size*2,UE_scheduling(u_).nCodewords),2);
%            end
           c(c<0) = 10^-2;
       end
   end
   
   methods (Abstract)
       % UE scheduling (to be implemented for each subclass
       UE_scheduling = scheduler_users(obj,subframe_corr,total_no_refsym,SyncUsedElements,UE_output,UE_specific_data)
   end
end 
