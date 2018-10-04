classdef UeSpecificEnodebData < handle
% This class stores resources that are used for each UE that is attached to
% this eNodeB. Eg. the HARQ processes.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       HARQ_processes
       current_HARQ_process
       free_HARQ_processes
   end

   methods
       function obj = UeSpecificEnodebData(maxStreams,HARQ_processes,max_rv_idx)
           obj.free_HARQ_processes = true(maxStreams,HARQ_processes);
           for cw_=1:maxStreams
               for harq_ = 1:HARQ_processes
                   obj.HARQ_processes{cw_,harq_} = network_elements.harqProcess(harq_,maxStreams,max_rv_idx);
               end
           end
       end
   end
end 