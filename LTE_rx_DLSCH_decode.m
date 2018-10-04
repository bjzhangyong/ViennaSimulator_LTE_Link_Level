function [a BER ACK ACK_b_binary C] = LTE_rx_DLSCH_decode(LTE_params,f,BS_signaling,UE,BS_ue_specific,genie_data,stream_index)
% Performs the decoding of a TB for the LTE LL simulator. Basically an implementation
% of TS 36.212. Naming convention also according to TS 36.212.
% [a BER ACK ] = LTE_rx_DLSCH_decode(f,BS_signaling,UE)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    f                   ... received LLRs (soft bits)
%           BS_signaling        ... struct containing parameters needed
%                                   for decoding (number of padding bits, 
%                                   TB size...)
%           UE                  ... struct that contains the number of 
%                                   turbo iterations that the user will
%                                   perform and some rate matching
%                                   params, such as the soft buffer size.
%
% output:   a                   ... decoded data bits
%           BER                 ... BER, as outputted by the turbo
%                                   decoder (using genie info). Actually
%                                   the number of erroneous bits per
%                                   turbo iteration
%           ACK                 ... ACK, calculated from the TB CRCs. A
%                                   true of false value stating if the TB
%                                   CRC check was correct (true) or
%                                   failed (false).
%           UE                  ... updated UE struct.
%
% Note: a BLER_b value containing the CRC checks of the individual code
% blocks is also computed, but neither used nor passed outside of the
% function.
%
% date of creation: 2008/08/11
% last changes: 2008/09/15 Bosanska     added input BS_ue_specific - [1 x 1]struct with tx user specific HARQ parameters
%               2010/11/03 JColom       removed LTE_rx_code_block_concatenation 

% Code block de-concatenation
e = mat2cell(f,BS_signaling.turbo_rate_matcher(stream_index).ek_sizes,1);

% Rate de-matching
d = LTE_rx_turbo_rate_matcher(LTE_params,e,BS_signaling,UE,stream_index);

% HARQ combining. Could have been embedded into the rate matcher, but it's less messy to put it separate
[ d_HARQ UE ] = LTE_rx_HARQ_combine(d,UE,BS_ue_specific,stream_index);

% Turbo decoding. Genie info is used to get BER values for each turbo iteration
[ c BER ] = LTE_rx_turbo_decode(LTE_params,d_HARQ,...
    UE.turbo_iterations,...
    BS_signaling.TB_segmentation(stream_index).C,...
    genie_data.bits_to_turboencode{stream_index});

% Code block desegmentation
[b ACK_b] = LTE_rx_code_block_desegmentation(c,BS_signaling.TB_segmentation(stream_index).F);

% Transport block CRC check
ACK = LTE_rx_check_crc(b,'24a');

% Get the Transport Block data
a = b(1:(end-24));

% Transcribe CB ACK to binary form
if isempty(ACK_b)
    ACK_b_binary = uint16(ACK);
else
    ACK_b_binary = sum(uint16(ACK_b).*2.^(uint16((1:length(ACK_b))-1)));
end

C = BS_signaling.TB_segmentation(stream_index).C;