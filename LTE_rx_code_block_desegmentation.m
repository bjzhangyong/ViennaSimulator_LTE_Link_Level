function [transport_block crc_check] = LTE_rx_code_block_desegmentation(code_blocks,F)
% Performs code block desegmentation and CRC checking according to TS 36.212, Section 5.1.2
% [transport_block crc_check] = LTE_rx_code_block_desegmentation(code_blocks,F)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   code_blocks         ... cell array containing the several code
%                                   blocks. code_block{1} is the first one
%           F                   ... number of filler bits used
% output:   transport_block     ... logical vector representing a transport
%                                   block
%           crc_check           ... for each code block, true or false,
%                                   depending on whether the crc was correct or
%                                   not
%
% date of creation: 2008/08/11
% last changes:

C = length(code_blocks);

if C==1
    % No extra CRC, return an empty segment CRC check
    transport_block = code_blocks{1}((F+1):end);
    crc_check = [];
else
    % First code block
    crc_check = false(1,C);
    crc_check(1) = LTE_rx_check_crc(code_blocks{1},'24b');
    transport_block = code_blocks{1}((F+1):(end-24)); % Take out filler bits
    % Rest of code blocks
    for i=2:C
        % Check CRCs
        crc_check(i) = LTE_rx_check_crc(code_blocks{i},'24b');
        transport_block = [transport_block code_blocks{i}(1:(end-24))];
    end
end