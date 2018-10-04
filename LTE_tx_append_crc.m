function data_and_crc = LTE_tx_append_crc(data,crc_type)
% Calculates and appends the indicated CRC to the input sequence data.
% [data_and_crc] = LTE_tx_append_crc(data,crc_type)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    data        ... logical vector representing the bits of the message to CRC
%           crc_type    ... which CRC to use. Can be {'16','24a','24b'}
% 
% output:   crc         ... calculated crc
%
% date of creation: 2008/08/11
% last changes:

%% Choose which CRC to use
switch crc_type
    case '16'
        crc_function = @LTE_common_crc16;
    case '24a'
        crc_function = @LTE_common_crc24a;
    case '24b'
        crc_function = @LTE_common_crc24b;
end

%% Calculate CRC. For CRC calculation it is assumed that filler bits, if present, have the value 0
data_bytes    = LTE_common_bit2byte(data);
crc_bytes     = crc_function(data_bytes);
crc           = LTE_common_byte2bit(crc_bytes);

%% Append CRC. According to TS 36.212, Section 5.1.1
data_and_crc = [ data crc ];
