classdef ueSpecificTraces < handle
% Class that stores all of the UE-scpecific traces
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

properties
    % Fields updated after every TTI
    ACK
    RBs_assigned
    rv_idx
    biterrors_coded
    biterrors_uncoded
    blocksize_coded
    blocksize_uncoded
    % CQI and RI mapping?
    FER_coded
    FER_uncoded
    throughput_coded
    throughput_uncoded
    throughput_useful
    freq_offset_est_error
    freq_offset_est_frac
    freq_offset_est_int
    freq_offset_est_res
    timing_offset_est
    
    % Aggregates calculated after the simulation is finisched
    BER_coded           % coded BER (each stream)
    BER_uncoded         % uncoded BER (each stream)
    BER_coded_overall   % coded BER (overall)
    BER_uncoded_overall % uncoded BER (overall)
    BLER                % BLER (each stream)
    BLER_overall        % BLER (overall)
    MSE_freq_offset     % MSE carrier frequency offset

    % Used to signal what entries are valid
    used_codewords
    
    % Codeblock traces
    used_codeblocks
    ACK_codeblocks
    C
    avg_CB_size
    
    % Log CQIs
    used_CQI
end

   methods
       % Class contructor. Data preallocation
       function obj = ueSpecificTraces(N_subframes,SNR_vector_length,maxStreams)
           obj.ACK                = false(N_subframes,SNR_vector_length,maxStreams);          % true/false -> ACK of the received subframes for BLER calculation
           obj.ACK_codeblocks     = zeros(N_subframes,SNR_vector_length,maxStreams,'uint16'); % stores the ACKs/NACKs of up to 16 codeblocks per stream. Stored in a binary form. i.e 1101100001000000
           obj.C                  = zeros(N_subframes,SNR_vector_length,maxStreams,'uint8');  % How many Codeblocks were used for each stream
           obj.avg_CB_size        = zeros(N_subframes,SNR_vector_length,maxStreams);          % Average Codeblock size for each codeword
           obj.RBs_assigned       = zeros(N_subframes,SNR_vector_length,'uint8');  % 0-255      -> number of assigned RBs in every TTI
           obj.rv_idx             = zeros(N_subframes,SNR_vector_length,maxStreams,'uint8');  % 0-255      -> redundancy version index of the received subframes
           obj.biterrors_coded    = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.biterrors_uncoded  = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.blocksize_coded    = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.blocksize_uncoded  = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.FER_coded          = false(N_subframes,SNR_vector_length,maxStreams);          % This is actually ~ACK
           obj.FER_uncoded        = false(N_subframes,SNR_vector_length,maxStreams);          %
           obj.throughput_coded   = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.throughput_uncoded = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.throughput_useful  = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % 0-4294967295
           obj.used_codeblocks    = zeros(N_subframes,SNR_vector_length,maxStreams,'uint32'); % How many codeblocks were used
           obj.used_codewords     = false(N_subframes,SNR_vector_length,maxStreams);          % What codewords were used
           obj.used_CQI           = zeros(N_subframes,SNR_vector_length,maxStreams,'uint8');
           obj.freq_offset_est_error = zeros(N_subframes,SNR_vector_length);
           obj.freq_offset_est_frac = zeros(N_subframes,SNR_vector_length);
           obj.freq_offset_est_int = zeros(N_subframes,SNR_vector_length);
           obj.freq_offset_est_res = zeros(N_subframes,SNR_vector_length);
           obj.timing_offset_est = zeros(N_subframes,SNR_vector_length);
       end
   end
end 
