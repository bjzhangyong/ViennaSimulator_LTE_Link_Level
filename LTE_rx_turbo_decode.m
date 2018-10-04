function [ output error_bits ] = LTE_rx_turbo_decode(LTE_params,LLR_input,iterations,C,varargin)
% LTE turbo decoder. Uses the Iterative Solutions Coded Modulation 
% Library SISO decoder (GPL).
% [ output error_bits ] = LTE_rx_turbo_decode(LLR_input,iterations,varargin)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   LLR_input           ... [1 x N]Cell array of doubles LLRs of the demapped 
%                                   received signal (soft bits).
%           iterations          ... [1 x 1]double Number of iterations for the turbo decoder
%           C                   ... [1 x 1]double Number of code blocks in
%                                   which this TB was divided. Needed to
%                                   correctly check the CRC.
%           data_bits(optional) ... [1 x N]logical If the original data bits
%                                   are inputted (genie info),
%                                   BER data will be outputted for each
%                                   turbo iteration. Cell array of bits to
%                                   encode,a s inputted to the encoder
% output:   output              ... [1 x N]Cell array of [1 x N]logicals.
%                                   Decoded data bits (hard bits)
%           error_bits          ... [1 x N]Cell array of [1 x iterations]logicals.
%                                   If the genie info is input, this
%                                   returns the number of erroneous bits
%                                   for each turbo iteration
%
% date of creation: 2008/08/11
% last changes:

length_LLR = length(LLR_input);
output = cell(length_LLR,1);
error_bits = cell(length_LLR,1);
for i=1:length_LLR
    if ~isempty(varargin)
        [ output{i} error_bits{i} ] = LTE_rx_turbo_decode_codeword(LTE_params,LLR_input{i},iterations,C,varargin{1}{i});
    else
        output{i} = LTE_rx_turbo_decode_codeword(LTE_params,LLR_input{i},iterations,C);
    end
end

function [ output error_bits ] = LTE_rx_turbo_decode_codeword(LTE_params,LLR_bits,iterations,C,varargin)
% LTE turbo decoder. Uses the Iterative Solutions Coded Modulation Library
% SISO decoder (GPL).
% (c) Josep Colom Ikuno, INTHFT
% josep.colom@nt.tuwien.ac.at
% www.nt.tuwien.ac.at
%
% input :   LLR_bits            ... LLRs of the demapped received signal (soft bits)
%           iterations          ... [1 x 1]double Number of iterations for
%                                   the turbo decoder
%           C                   ... [1 x 1]double Number of code blocks in
%                                   which this TB was divided. Needed to
%                                   correctly check the CRC.
%           data_bits(optional) ... [1 x N]logical If the original data bits
%                                   are inputted (genie info),
%                                   BER data will be outputted for each
%                                   turbo iteration. Cell array of bits to
%                                   encode,a s inputted to the encoder
% output:   output              ... [1 x N]logical Decoded data bits (hard bits)
%           error_bits          ... [1 x iterations]logical If the genie info is input, this
%                                   returns the number of erroneous bits
%                                   for each turbo iteration

% Generator polynomial
g_turbo = [1 0 1 1
    1 1 0 1];

% Decoder type (code tested only with the max-log-map configuration, not with any other decoder type)
decoder_type = LTE_params.UE_config.decoder_type_idx;

% decoder_type = 0; % For linear approximation to log-MAP (DEFAULT)
% decoder_type = 1; % For max-log-MAP algorithm (i.e. max*(x,y) = max(x,y) )
% decoder_type = 2; % For Constant-log-MAP algorithm
% decoder_type = 3; % For log-MAP, correction factor from small nonuniform table and interpolation
% decoder_type = 4; % For log-MAP, correction factor uses C function calls (slow)

nsc_flag = 0; % Recursive-systematic code
% Get Code Block size
code_block_size = (size(LLR_bits,1)*size(LLR_bits,2)-12)/3;
% Get interleaver mapping
interleaver_mapping = LTE_common_turbo_encoder_generate_interleaving_mapping(LTE_params,code_block_size);

% Add tail bits
y0_1 = [ LLR_bits(1,1:end-4) LLR_bits(1,end-3) LLR_bits(3,end-3) LLR_bits(2,end-2) ];
y1   = [ LLR_bits(2,1:end-4) LLR_bits(2,end-3) LLR_bits(1,end-2) LLR_bits(2,end-2) ];

% y0_2 = [ LTE_common_soft_bit_interleaver(LLR_bits(1,1:end-4),interleaver_mapping,1) LLR_bits(1,end-1) LLR_bits(3,end-1) LLR_bits(2,end)   ];
y0_2 = zeros(size(y1)); % Initialize systematic part of second decoder with zero conforming CML. Thanks for finding this error go to Klaus Hueske :)
y2   = [ LLR_bits(3,1:end-4) LLR_bits(2,end-1) LLR_bits(1,end)   LLR_bits(3,end)   ];

SISO1_input = reshape([y0_1;y1],1,length(y0_1)+length(y1));
SISO2_input = reshape([y0_2;y2],1,length(y0_2)+length(y2));

%% Initialise APP from SISO1
SISO1_APP = zeros(1,code_block_size);

%% If there is genie info, check the BER of the uncoded (systematic) bits
if length(varargin)>=1
    error_bits = zeros(1,iterations+1);
    received_bits = LTE_rx_hard_decision(LLR_bits(1,1:end-4));
    error_bits(1)=sum(xor(received_bits,varargin{1,1}));
else
    error_bits = [];
end

%% Choose which CRC to use to check the code blocks
% If the TB is just one code block, then CRC24a is used. If not, CRC24b.
switch C
    case 1
        crc_type = '24a';
    otherwise
        crc_type = '24b';
end

%% Start iterations
for turbo_iteration = 1:iterations
    % First decoder
    SISO1_output = LTE_rx_siso_decode(SISO1_APP,SISO1_input,g_turbo,nsc_flag,decoder_type);
    % Second decoder
    SISO2_APP = LTE_common_soft_bit_interleaver(SISO1_output-SISO1_APP,interleaver_mapping,1);
    SISO2_output = LTE_rx_siso_decode(SISO2_APP,SISO2_input,g_turbo,nsc_flag,decoder_type);
    SISO1_APP = LTE_common_soft_bit_interleaver(SISO2_output-SISO2_APP,interleaver_mapping,0);
    
    % Get decoded bits so far
    decoded_bits = LTE_common_bit_interleaver(LTE_rx_hard_decision(SISO2_output),interleaver_mapping,0);
    
    % Compute BER, if genie bits are present
    if length(varargin)>=1
        error_bits(turbo_iteration+1)=sum(xor(decoded_bits,varargin{1,1}));
    end
    
    % If CRC is correct, halt iterations
    if LTE_rx_check_crc(decoded_bits,crc_type)
        break;
    end
end

%% De-interleave the output of the second decoder
output = LTE_common_bit_interleaver(LTE_rx_hard_decision(SISO2_output),interleaver_mapping,0);