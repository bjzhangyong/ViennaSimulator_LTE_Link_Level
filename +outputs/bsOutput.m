classdef bsOutput < handle
% Output of the eNodeB. Contains data, signaling information and genie
% data also.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       y_tx                     % The sent data
       UE_signaling             % Signaling for each UE
       genie                    % Genie information
       cell_genie               % Genie information related to the whole cell (eg. the Subcarrier SNR)

       
   end

   methods
       function obj = bsOutput(nUE,nRB,maxStreams)
           obj.UE_signaling = outputs.ueSignaling;
           obj.genie = outputs.genieInformation;
           for u_=1:nUE
               obj.UE_signaling(u_) = outputs.ueSignaling;
               obj.genie(u_) = outputs.genieInformation;
           end
       end
       
       function clear(obj)
           obj.y_tx = [];
           for u_=1:length(obj.UE_signaling)
               obj.UE_signaling(u_).clear;
               obj.genie(u_).clear;
           end

       end
   end
end 
