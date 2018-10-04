classdef UE_scheduling_params < handle
% Scheduling parameters for a specific UE
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
	   tx_mode             % Assigned TX mode:
	                       %   1: Single Antenna, 2: Transmit Diversity,
	                       %   3: Open Loop Spatial Multiplexing, 4: Closed Loop SM
						   
       nLayers             % Number of layers used during this TTI for this UE.
	                       % With this value and the transmission mode it is
						   % possible to get how many codewords there are. Anyway,
						   % for receiver simplicity, I will also signal it.
						   
	   nCodewords          % Number of codewords used this TTI
       
       cqi                 % One CQI per codeword
       UE_mapping          % RBs assigned to this UE. Accessed as UE_mapping(freq_index,slot_index)
       assigned_RBs        % How many RBs where allocated for each stream for this UE
       CQI_params          % Parameters related to the CQI used
       N_coded_bits        % Number of coded bits
       N_data_bits         % Number of data bits (obtained by applying the Effective Code Rate).
       N_used_bits         % Number of actually used bits - there might be too little data to fill a full RB
       rv_idx              % HARQ Redundancy version index to use
       HARQ_process_id     % ID of the current HARQ process for each stream
       
       % MIMO parameters
       PMI
       codebook_index
       CDD
       PRE
       freq_indices
       slot_indices
       V
       U       
   end

   methods
%        function obj = UE_scheduling_params(cqi_params)
%            obj.CQI_params = struct(cqi_params);
%        end
   end
end 
