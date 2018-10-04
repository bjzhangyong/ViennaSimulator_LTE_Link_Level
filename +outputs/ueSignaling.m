classdef ueSignaling < handle
% DL-SCH signaling for each UE.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       % Channel coding
       turbo_rate_matcher  % Signaling related to turbo rate matching
       TB_size             % TB size
       TB_segmentation     % Information about TB segmentation
       turbo_encoder       % Turbo encoder related info
       
       % Scheduling
	   MCS_and_scheduling % Information related to modulation and coding (CQI-related), scheduling and precoding
   end

   methods
       function clear(obj)
           obj.turbo_rate_matcher = [];
           obj.TB_size = [];
           obj.TB_segmentation = [];
           obj.turbo_encoder = [];
           obj.MCS_and_scheduling = [];
       end
   end
end 
