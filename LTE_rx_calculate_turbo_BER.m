function BER = LTE_rx_calculate_turbo_BER(BER,BS_signaling,BER_error_bits_turbo)
% Processes the data from the current iteration and adds it to the BER
% struct containing the accumulated BER data so far.
% [BER] = LTE_rx_calculate_turbo_BER(BER,BS_signaling,BER_error_bits_turbo)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    BER                   ... struct containing the accumulated BER data so far
%           BS_signaling          ... struct containing the signaling from the BS
%           BER_error_bits_turbo  ... Erroneous bits in the current TB
%
% output:   BER                   ... updated BER struct
%
% date of creation: 2008/08/11
% last changes:

global LTE_params;

total_bits = sum(BS_signaling.turbo_encoder.encoded_CB_sizes);
total_erroneous_bits = zeros(1,LTE_params.turbo_iterations+1);

for i=1:length(BER_error_bits_turbo)
    total_erroneous_bits = total_erroneous_bits + BER_error_bits_turbo{i};
end

BER.total_bits = BER.total_bits + total_bits;
BER.erroneous_bits = BER.erroneous_bits + total_erroneous_bits;