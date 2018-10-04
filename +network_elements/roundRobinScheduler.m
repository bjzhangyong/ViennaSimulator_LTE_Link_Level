classdef roundRobinScheduler < network_elements.lteScheduler
% A round robin scheduler (equally schedule every user). Please note the
% following:
%  - For a static scheduler: since all the scheduling info is generated at
%    object creation, the UE allocation will be time invariant. Thus, if
%    the number of RBs is not an integer multiple of the number of UEs,
%    some RBs will be left unused.
%  - For a dynamic scheduler: still to be implemented :P
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
   end

   methods
       function obj = roundRobinScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params)
           
           % Fill in basic parameters (handled by the superclass constructor)
           obj = obj@network_elements.lteScheduler(RB_grid_size,Ns_RB,UEs_to_be_scheduled,scheduler_params,CQI_params);

           switch scheduler_params.assignment
               case 'static'
                   obj.static_scheduler = true;
                   number_of_RBs_per_UE = floor(RB_grid_size*2 / obj.nUEs);
                   
                   % Get a vector of scheduling params (one for each UE)
                   % initialized to the values that we want
                   obj.UE_static_params = obj.get_initialized_UE_params(scheduler_params,CQI_params);
                   
                   % Fill in the RB allocation grid for each user (and codeword)
                   %UE_mapping_all_UEs = zeros(RB_grid_size,2,obj.maxCodewords);
                   UE_mapping_all_UEs = zeros(RB_grid_size,2);
                   for u_=1:obj.nUEs
                       % NOTE: Same RB assignment for both codewords.
                       UE_RBs = 1 + number_of_RBs_per_UE*(u_-1) : number_of_RBs_per_UE + number_of_RBs_per_UE*(u_-1);
                       cw_RB_grid = UE_mapping_all_UEs;
                       cw_RB_grid(UE_RBs) = u_;
                       UE_mapping_all_UEs = cw_RB_grid;
%                        for cw_=2:scheduler_params.nCodewords(u_)
%                            UE_mapping_all_UEs(:,cw_) = UE_mapping_all_UEs(:,1);
%                        end
                   end

                   % Assign the static scheduling parameters for each user
                   for u_=1:obj.nUEs
                       obj.UE_static_params(u_).UE_mapping = (UE_mapping_all_UEs==u_);
                       obj.UE_static_params(u_).assigned_RBs = squeeze(sum(sum(obj.UE_static_params(u_).UE_mapping,1),2));
                   end
               case 'dynamic'
                   obj.static_scheduler = false;
                   % NOTE: one way of implementing a dynamic scheduler
                   % could be to actually have the "static" scheduler but
                   % with some more degrees of freedom. Eg. Fix transmit
                   % mode but allow for CQI variations. Implementing all of
                   % the degrees of freedom at the same time could lead to
                   % too many parameters to keep track of for an earlier
                   % scheduler algorithm.
                   error('dynamic scheduler not yet implemented!!!!');
               otherwise
                   error('only static or dynamic schedulers supported');
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
           if obj.static_scheduler
               UE_scheduling = obj.copy(obj.UE_static_params);
%                UE_scheduling = obj.UE_static_params;
           else   
               error('Dynamic scheduler not yet implemented.');
           end
           obj.calculate_allocated_bits(UE_scheduling,subframe_corr,total_no_refsym,SyncUsedElements,PBCHsyms);
       end
       
        function new = copy(obj,in)
           % Instantiate new object of the same class.
           for u_=1:obj.nUEs
               new(u_) = feval(class(in));
               
               % Copy all non-hidden properties.
               p = properties(in);
               for i = 1:length(p)
                   new(u_).(p{i}) = in(u_).(p{i});
               end
           end
       end
   end
end 
