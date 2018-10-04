classdef eNodeB < handle
% Class that represents an eNodeB.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       NIDcell              % Cell identity
       UE_specific          % UE specific data, such as the HARQ TX buffer
       user_count           % How many users there are
       nTX                  % Number of transmit antennas
       scheduler            % Where we will store this eNodeB's scheduler
       AtPort               % supported antenna ports 
       nAtPort              % number of antenna ports
       clock                % So the eNodeB is aware in which TTI he is
       HARQ_process_count   % Number of HARQ processes
   end

   methods
       % Class constructor
       function obj = eNodeB(NID,nTx,nUE,AtPort,nAtPort,maxStreams,HARQ_processes,max_rv_idx)
           obj.NIDcell = NID;
           obj.UE_specific = network_elements.UeSpecificEnodebData(maxStreams,HARQ_processes,max_rv_idx);
           for uu=1:nUE
               obj.UE_specific(uu) = network_elements.UeSpecificEnodebData(maxStreams,HARQ_processes,max_rv_idx);
           end
           obj.nTX = nTx;
           obj.user_count = nUE;
           obj.AtPort = AtPort;
           obj.nAtPort = nAtPort;
           obj.HARQ_process_count = HARQ_processes;
       end
       
       % Advance the HARQ process that is being used
       function update_current_HARQ_process(obj,UE_output)

           for u_=1:obj.user_count
               
               max_codewords = size(obj.UE_specific(u_).HARQ_processes,1);
               % empty HARQ process means one of the first ones, where no ACKs are present.
               if isempty(UE_output(u_).HARQ_process)
                   % Set the current HARQ process as the first free one
                   for cw_=1:max_codewords
                       if obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs > 0
                           HARQ_idx(cw_) = find(obj.UE_specific(u_).free_HARQ_processes(cw_,:),1,'first');
                           obj.UE_specific(u_).free_HARQ_processes(cw_,HARQ_idx(cw_)) = false;
                           obj.UE_specific(u_).current_HARQ_process(cw_) = obj.UE_specific(u_).HARQ_processes{cw_,HARQ_idx(cw_)};
                           obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = 0;
                       else
                           HARQ_idx(cw_) = obj.UE_specific(u_).current_HARQ_process(cw_).id;
                           obj.UE_specific(u_).free_HARQ_processes(cw_,HARQ_idx(cw_)) = false;
                       end
                   end
               else
                   % Mark HARQ process from which the current ACK has been received as free
                   UE_scheduled      = UE_output(u_).UE_scheduled;
                   feedback_HARQ_idx = [ UE_output(u_).HARQ_process(:);zeros(max_codewords-length(UE_output(u_).HARQ_process(:)),1)]; % Processes from which the current ACK comes
                   
                   for cw_=1:length(UE_scheduled)
                       
                       % Mark newly-free HARQ processes
                       if UE_scheduled(cw_)
                           if feedback_HARQ_idx(cw_)~=0
                               obj.UE_specific(u_).free_HARQ_processes(cw_,feedback_HARQ_idx(cw_)) = true;
                           end
                       end
                       
                       % TX HARQ process
                       % If the current HARQ process was used (ie. scheduled), move the current HARQ process forward and update rv_idx
                       % Otherwise do nothing
                       if obj.UE_specific(u_).current_HARQ_process(cw_).assigned_RBs > 0
                           next_free_HARQ_idx(cw_) = find(obj.UE_specific(u_).free_HARQ_processes(cw_,:),1,'first');
                           obj.UE_specific(u_).current_HARQ_process(cw_) = obj.UE_specific(u_).HARQ_processes{cw_,next_free_HARQ_idx(cw_)};
                           current_rv_idx = obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx;
                           max_rv_idx     = obj.UE_specific(u_).current_HARQ_process(cw_).max_rv_idx;
                           if isempty(current_rv_idx)
                               obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = 0;
                           else
                               if UE_output(u_).ACK(cw_)
                                   obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = 0;
                                   if strcmp(obj.scheduler.UEs(u_).traffic_model.type,'voip') || strcmp(obj.scheduler.UEs(u_).traffic_model.type,'video') || strcmp(obj.scheduler.UEs(u_).traffic_model.type,'gaming')
                                        for pp = 1:length(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts) % acknowledge all packet parts and remove them from the buffer
%                                             [packet_done,packet_id] = obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).acknowledge_packet_part(true);
%                                             if packet_done
%                                                 obj.scheduler.UEs(u_).traffic_model.remove_packet(packet_id,true);
%                                             end
                                            if obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).data_packet_id
                                                packet_ind = obj.scheduler.UEs(u_).traffic_model.get_packet_ids == obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).data_packet_id;
                                                if sum(packet_ind)
                                                    [packet_done,packet_id] = obj.scheduler.UEs(u_).traffic_model.packet_buffer(packet_ind).acknowledge_packet_part(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).id,true);
    %                                                 [packet_done,packet_id] = obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).acknowledge_packet_part(true);
                                                    if packet_done && packet_id
                                                        obj.scheduler.UEs(u_).traffic_model.remove_packet(packet_id,true);
                                                    end
                                                end
                                            end
                                        end
                                   end
                               else
                                   obj.UE_specific(u_).current_HARQ_process(cw_).rv_idx = mod(current_rv_idx+1,max_rv_idx+1);
                                   if strcmp(obj.scheduler.UEs(u_).traffic_model.type,'voip') || strcmp(obj.scheduler.UEs(u_).traffic_model.type,'video') || strcmp(obj.scheduler.UEs(u_).traffic_model.type,'gaming')
                                        if current_rv_idx == max_rv_idx % if maximum of retransmissions is obtained --> delete the packet parts
                                            for pp = 1:length(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts)
%                                                 [packet_done,packet_id] = obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).acknowledge_packet_part(false);
%                                                 if packet_done
%                                                     obj.scheduler.UEs(u_).traffic_model.remove_packet(packet_id,false);
%                                                 end
                                                if obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).data_packet_id
                                                    packet_ind = obj.scheduler.UEs(u_).traffic_model.get_packet_ids == obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).data_packet_id;
                                                    if sum(packet_ind)
                                                        [packet_done,packet_id] = obj.scheduler.UEs(u_).traffic_model.packet_buffer(packet_ind).acknowledge_packet_part(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).id,false);
                                                        if packet_done && packet_id
                                                            obj.scheduler.UEs(u_).traffic_model.remove_packet(packet_id,false);
                                                        end
                                                    end
                                                end
                                            end
                                        else % else restore them for retransmission
                                            for pp = 1:length(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts)
%                                                 obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).restore_packet_part;
%                                                 delete(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp));
                                                if obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).data_packet_id
%                                                     obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).restore_packet_part;
%                                                     delete(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp));
                                                    packet_ind = obj.scheduler.UEs(u_).traffic_model.get_packet_ids == obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).data_packet_id;
                                                    [packet_done,packet_id] = obj.scheduler.UEs(u_).traffic_model.packet_buffer(packet_ind).restore_packet_part(obj.UE_specific(u_).current_HARQ_process(cw_).packet_parts(pp).id);
%                                                     if packet_done && packet_id
%                                                         obj.traffic_model.remove_packet(packet_id,true);
%                                                     end
                                                end
                                            end
                                        end
                                   end
                               end
                           end
                           obj.UE_specific(u_).free_HARQ_processes(cw_,next_free_HARQ_idx(cw_)) = false;
                       else
                           % Current HARQ process is still the same
                       end
                   end
               end
           end
       end
       
       % Reset the HARQ process index
       function reset_HARQ_process_index(obj)
           for u_=1:obj.user_count
               % Set the HARQ processes for all the codewords
               obj.UE_specific(u_).current_HARQ_process = [obj.UE_specific(u_).HARQ_processes{:,1}];
               % Reset the rest of the HARQ processes
               for cw_=1:size(obj.UE_specific(u_).HARQ_processes,1)
                   for harq_=1:size(obj.UE_specific(u_).HARQ_processes,2)
                       obj.UE_specific(u_).HARQ_processes{cw_,harq_}.reset;
                       obj.UE_specific(u_).free_HARQ_processes(cw_,:) = ones(size(obj.UE_specific(u_).free_HARQ_processes(cw_,:)));
                   end
               end
           end
       end
   end
end 
