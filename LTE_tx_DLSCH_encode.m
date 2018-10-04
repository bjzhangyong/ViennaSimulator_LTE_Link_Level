function f = LTE_tx_DLSCH_encode(LTE_params,a,UE_signaling,UE_genie,UE,stream_index)
% Performs the coding of a TB for the LTE LL simulator. Basically an implementation
% of TS 36.212. Naming convention also according to TS 36.212.
% [f coding_params ] = LTE_tx_DLSCH_encode(a,scheduler_params,UE)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    a                     ... data bits
%           UE_signaling          ... Contains the object that will be sent as signaling to the UE.
%           UE                    ... UE information. We just need N_IR
%                                     params, such as the soft buffer size.
%           UE_genie              ... Object where the genie information will be written
%           UE                    ... Object representing the UE
%           stream_index          ... Codeword index
%
% output:   f                     ... coded bits
%
% date of creation: 2008/08/11
% last changes: 2008/09/15  Bosanska    changed input [1x1]struct ue_specific_params according to changed BS_output
%                                       added input [1x1]struct scheduler_params according to changed BS_output
%                                       added outputs [1x1]struct
%                                       coding_params and genie according
%                                       to changed BS_output

% Store the original data sequence, for genie-driven calculations, such as BER.
UE_genie.data_bits{stream_index} = a;

% Transport Block CRC attachment
b = LTE_tx_append_crc(a,'24a');
UE_signaling.TB_size(stream_index) = length(b); % Store the TB size (signal it to the receiver)

% Code Block segmentation and CRC attachment
c = LTE_tx_code_block_segmentation(LTE_params,b,UE_signaling,stream_index);

% Channel Coding (Turbo encoder)
d = LTE_tx_turbo_encode(LTE_params,c,UE_signaling,stream_index);
UE_genie.bits_to_turboencode{stream_index} = c;

% Rate matching
e = LTE_tx_turbo_rate_matcher(LTE_params,d,UE_signaling,UE,stream_index);

% Code block concatenation
f = LTE_tx_code_block_concatenation(e);
UE_genie.sent_bits{stream_index} = f;
