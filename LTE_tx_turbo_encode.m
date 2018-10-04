function coded_blocks = LTE_tx_turbo_encode(LTE_params,code_blocks,UE_signaling,stream_index)
% LTE turbo encoder.
% [coded_blocks BS_signaling] = LTE_tx_turbo_encode(code_blocks,BS_signaling)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    code_blocks   ... cell array (vector) containing N code blocks
%                             (logical vectors).
%           BS_signaling  ... where some signaling information will be
%                             stored.
% 
% output:   coded_blocks  ... coded blocks. Cell array of encoded code 
%                             blocks. Each coded block, since this is 
%                             a 1/3 systematic code, has 3 rows, one 
%                             ouputted coded bit.
%           BS_signaling  ... in order to get the data back, BS_signaling
%                             must be returned.
%
% Each constituent convolutional encoder uses a generator poly: g0/g1,
% where (TS 36.212, Section 5.1.3.2.1):
%   g0=1+D^2+D^3
%   g1=1+D+D^3
%
% date of creation: 2008/08/11
% last changes:

C = UE_signaling.TB_segmentation(stream_index).C;
coded_blocks = cell(1,C);
for i=1:C
    coded_blocks{i} = LTE_tx_turbo_encode_code_block(LTE_params,code_blocks{i},i,UE_signaling.TB_segmentation(stream_index));
    UE_signaling.turbo_encoder(stream_index).encoded_CB_sizes(i) = length(coded_blocks{i});
end

function coded_bits = LTE_tx_turbo_encode_code_block(LTE_params,input_bits,cb_idx,TB_segmentation_info)
% LTE turbo encoder.
% (c) Josep Colom Ikuno, INTHFT
% josep.colom@nt.tuwien.ac.at
% www.nt.tuwien.ac.at
%
% input:    input_bits ... logical vector representing the input bits. This
%                          is a CODEWORD, already after segmentation.
% 
% output:   coded_bits ... coded bits. This is a 1/3 systematic code, so
%                           the output has 3 row, one for each output bit.
%
% Each constituent convolutional encoder uses a generator poly: g0/g1,
% where (TS 36.212, Section 5.1.3.2.1):
%   g0=1+D^2+D^3
%   g1=1+D+D^3

% Get Code Block size
code_block_size = length(input_bits);

% Get interleaver mapping
interleaver_mapping = LTE_common_turbo_encoder_generate_interleaving_mapping(LTE_params,code_block_size);

% Interleave one input
input_to_e1 = input_bits;
input_to_e2 = LTE_common_bit_interleaver(input_bits,interleaver_mapping,1);

% The same generator poly for the 2 Conv constituent encoders
generator_poly = [
    1 0 1 1 
    1 1 0 1 ];

% We are using a RSC code, this flag signals it
nsc_flag = 0;

% Encode the input sequence. The output already contains the terminating bits
encoded_bits_1 = LTE_tx_convolutional_encoder(input_to_e1,generator_poly,nsc_flag);
encoded_bits_2 = LTE_tx_convolutional_encoder(input_to_e2,generator_poly,nsc_flag);

data_size = length(input_bits);
tail_size = size(encoded_bits_1,2)-data_size;

% Add the encoded data bits
encoded_bits = [
    encoded_bits_1(1,1:data_size) % Systematic bits:       d_k^0
    encoded_bits_1(2,1:data_size) % Parity from encoder 1: d_k^1
    encoded_bits_2(2,1:data_size) % Parity from encoder 2: d_k^2
    ];

% If the code block to be encoded is the 0-th code block and the number of filler bits is greater than zero,
% i.e., F > 0, then the encoder shall set c_k=0, k=0,…,(F-1) at its input and shall set  d^0_k=<NULL>, 0,…,(F-1) and
% d^1_k=<NULL>, 0,…,(F-1) at its output.
if cb_idx==1
    if TB_segmentation_info.F>0
        encoded_bits([1 2],1:TB_segmentation_info.F) = LTE_params.null_bit;
    end
end

% Add the tailing bits according to what TS36.212, Section 5.1.3.2.2 says
tail_bits = [
    encoded_bits_1(1,data_size+1) encoded_bits_1(2,data_size+2) encoded_bits_2(1,data_size+1) encoded_bits_2(2,data_size+2)
    encoded_bits_1(2,data_size+1) encoded_bits_1(1,data_size+3) encoded_bits_2(2,data_size+1) encoded_bits_2(1,data_size+3)
    encoded_bits_1(1,data_size+2) encoded_bits_1(2,data_size+3) encoded_bits_2(1,data_size+2) encoded_bits_2(2,data_size+3)
    ];

coded_bits = [ encoded_bits tail_bits ];
