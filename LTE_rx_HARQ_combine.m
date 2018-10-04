function [ d_HARQ UE ] = LTE_rx_HARQ_combine(d,UE,BS_ue_specific,stream_index)
% Combines the demapped and rate de-matched LLR soft bitstream with the
% previous tranmissions, stored in a buffer in UE. Updates the buffer also.
% [ d_HARQ UE ] LTE_rx_HARQ_combine(d,UE)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) by INTHFT
% www.nt.tuwien.ac.at
%
% input :   d                   ... [3 x N]double  - output of the rate de-matcher
%           UE                  ... [1 x 1]struct  - User equipments capabilities from LTE_load_parameters.m
% output:   d_HARQ              ... [3 x N]double  - output of the rate, de-matcher with the additional info that the HARQ can provide
%           UE                  ... [1 x 1]struct  - User equipments capabilities from LTE_load_parameters.m,
%                                                    updates the receive soft buffer
%
% date of creation: 2008/08/21
% last changes: 2008/09/15  Bosanska     added input BS_ue_specific - [1 x 1]struct with tx user specific HARQ parameters

rv_idx = BS_ue_specific.current_HARQ_process(stream_index).rv_idx;
HARQ_process_idx = BS_ue_specific.current_HARQ_process(stream_index).id;

% If this is the first tx, then fill in the soft rx-buffer
if rv_idx==0
    UE.HARQ_rx_soft_buffer{HARQ_process_idx,stream_index} = d;
    d_HARQ = d;
else
    % If not, then add to the buffer each codewords (max. 6144 bits)
    for ii=1:length(d)
        UE.HARQ_rx_soft_buffer{HARQ_process_idx,stream_index}{ii} = UE.HARQ_rx_soft_buffer{HARQ_process_idx,stream_index}{ii} + d{ii};
    end
    d_HARQ = UE.HARQ_rx_soft_buffer{HARQ_process_idx,stream_index};
end

