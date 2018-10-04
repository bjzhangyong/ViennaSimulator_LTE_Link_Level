% Compile LTE Link Level Simulator MEX files
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2008/08/11
% last changes:

DEBUG = false;
clear mex;
%% MEX files that have a debug version
if DEBUG
    fprintf('MEX-ing files WITH debug output\n');
    mex ./C-source/crc16.c      -DDEBUG -output LTE_common_crc16
    mex ./C-source/crc24a.c     -DDEBUG -output LTE_common_crc24a
    mex ./C-source/crc24b.c     -DDEBUG -output LTE_common_crc24b
    mex ./C-source/ConvEncode.c -DDEBUG -output LTE_tx_convolutional_encoder
    mex ./C-source/byte2bit.c   -DDEBUG -output LTE_common_byte2bit
    mex ./C-source/bit2byte.c   -DDEBUG -output LTE_common_bit2byte
else
    fprintf('MEX-ing files WITHOUT debug output\n');
    mex ./C-source/crc16.c      -output LTE_common_crc16
    mex ./C-source/crc24a.c     -output LTE_common_crc24a
    mex ./C-source/crc24b.c     -output LTE_common_crc24b
    mex ./C-source/ConvEncode.c -output LTE_tx_convolutional_encoder
    mex ./C-source/byte2bit.c   -output LTE_common_byte2bit
    mex ./C-source/bit2byte.c   -output LTE_common_bit2byte
end

%% MEX files without a debug version
mex ./C-source/Bit_Interleaver.c                                      -output LTE_common_bit_interleaver
mex ./C-source/Soft_Bit_Interleaver.c                                 -output LTE_common_soft_bit_interleaver
mex ./C-source/SisoDecode.c                                           -output LTE_rx_siso_decode
mex ./C-source/Hard_decision.c                                        -output LTE_rx_hard_decision
%mex ./C-source/Circular_Buffer_Mapping.c                              -output LTE_common_turbo_rate_matching_circular_buffer_mapping
mex ./C-source/Turbo_rate_matching_bit_selection.c                    -output LTE_common_turbo_rate_matching_bit_selection_and_pruning
%mex ./C-source/Turbo_rate_matcher_bit_selection_and_pruning_mapping.c -output LTE_common_turbo_rate_matcher_bit_selection_and_pruning_mapping
mex ./C-source/soft_sd2.cpp                                           -output LTE_rx_soft_sd2
mex ./C-source/Gold_sequence_generation.c                             -output LTE_common_gen_gold_sequence