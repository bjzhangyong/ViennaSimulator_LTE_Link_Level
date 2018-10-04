classdef uplinkChannel < handle
% Class that represents the uplink channel. Basically a buffer containing N
% feefdback outputs, which delays the UE feedback a specified amount of
% TTIs.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       insert_idx
       extract_idx
       buffer_length
       buffer
   end

   methods
       % Class constructor. Delay >= 1 (minimum possible delay is 1 TTI)
       function obj = uplinkChannel(delay,nUEs)
           obj.insert_idx  = delay;
           obj.extract_idx = 1;
           
           obj.buffer_length = delay;
           if delay == 0
               obj.buffer_length = 1;
           end
           % Fill in buffer with the correct object type
           obj.buffer = outputs.ueOutput;

           for u_ = 1:nUEs             % Over UEs
               for d_ = 1:obj.buffer_length    % Over delay
                   obj.buffer(d_,u_) = outputs.ueOutput;
               end
           end 
       end
       
       % Insert UE feedbacks
       function insert_feedback(obj,UE_feedbacks)
           obj.buffer(obj.insert_idx,:) = UE_feedbacks;
       end
       
       % Extract UE feedbacks and advances the insert and extraction
       % indices one position
       function UE_feedbacks = receive_feedback(obj)
           UE_feedbacks = obj.buffer(obj.extract_idx,:);
           obj.advance_idxs;
       end
       
       % Advance the insert and extract indices one position
       function advance_idxs(obj)
           obj.insert_idx  = obj.insert_idx  + 1;
           obj.extract_idx = obj.extract_idx + 1;
           
           if obj.insert_idx > obj.buffer_length
               obj.insert_idx = 1;
           end
           
           if obj.extract_idx > obj.buffer_length
               obj.extract_idx = 1;
           end
       end
   end
end 
