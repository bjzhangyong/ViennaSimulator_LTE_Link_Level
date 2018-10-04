classdef ueOutput < handle
% UE output, including ACK, subframe rv_idx, CQI...
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       ACK
       ACK_codeblocks
       C
       UE_scheduled
       rv_idx
       CQI_feedback
       rx_data_bits
       rx_coded_bits
       RI
       PMI
       CQI
       CQI_bar
       HARQ_process
       channel_estimation_error
       freq_offset_est
       timing_offset_est
       PE_signal_power_subframe % Post Equalization Signal Power Matrix
       PE_noise_power_subframe  % Post Equalization Noise Power Matrix
       Signal_plus_noise_power
       Noise_power
   end

   methods
       function clear(obj)
           obj.ACK = [];
           obj.ACK_codeblocks = [];
           obj.C = [];
           obj.UE_scheduled = [];
           obj.rv_idx = [];
           obj.CQI_feedback = [];
           obj.rx_data_bits = [];
           obj.rx_coded_bits = [];
           obj.RI = [];
           obj.PMI = [];
           obj.HARQ_process = [];
           obj.channel_estimation_error = [];
           obj.freq_offset_est = [];
           obj.timing_offset_est = [];
           obj.CQI = [];
           obj.CQI_bar = 0;
           obj.PE_signal_power_subframe = [];
           obj.PE_noise_power_subframe = [];
           obj.Signal_plus_noise_power = [];
           obj.Noise_power = [];
       end
   end
end 
