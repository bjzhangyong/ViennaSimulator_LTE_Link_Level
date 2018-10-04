function crc_correct = LTE_rx_check_crc(data,crc_type)
% Checks whether the CRC is correct for the given data sequence.
% [crc_correct] = LTE_rx_check_crc(data,crc_type)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    data        ... [1 x N]logical The bits of the
%                           message you want to check (msg+crc).
%           crc_type    ... [string] Which CRC to use. Can be {'16','24a','24b'}
% 
% output:   crc_correct ... [1 x logical] Whether the CRC was correct.
%
% date of creation: 2008/08/11
% last changes:

data_length = length(data);

%% Choose which CRC to use and segment data
switch crc_type
    case '16'
        crc_function = @LTE_common_crc16;
        message_data = data(1:data_length-16);
        message_crc = data(data_length-15:end);
    case '24a'
        crc_function = @LTE_common_crc24a;
        message_data = data(1:data_length-24);
        message_crc = data(data_length-23:end);
    case '24b'
        crc_function = @LTE_common_crc24b;
        message_data = data(1:data_length-24);
        message_crc = data(data_length-23:end);
end

%% Check CRC
data_bytes = LTE_common_bit2byte(message_data);
message_crc_bytes = LTE_common_bit2byte(message_crc);

calculated_crc_bytes = crc_function(data_bytes);

switch min(message_crc_bytes==calculated_crc_bytes)
    case 0
        crc_correct = false;
    otherwise
        crc_correct = true;
end