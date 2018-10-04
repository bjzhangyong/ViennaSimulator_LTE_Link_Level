function [bin_BLER bins_begin] = LTE_common_compute_BLER(EESM_vector,ACK_vector,num_bins)
% From a list of SINRs and ACKs, it bins the results and returns back the
% corresponding BLER curve.
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

% Ignore NaNs and Infs
finites_idx = isfinite(EESM_vector);
EESM_vector = EESM_vector(finites_idx);
ACK_vector = ACK_vector(finites_idx);

min_SINR = min(EESM_vector);
max_SINR = max(EESM_vector);
bin_size = (max_SINR-min_SINR)/num_bins;

bins_begin = min_SINR:bin_size:max_SINR-bin_size;
%bin_middles = bins_begin + bin_size/2;
bin_sizes  = zeros(num_bins,1);
bin_BLER   = nan(num_bins,1);

for bin_idx = 1:num_bins
    init_SINR = bins_begin(bin_idx);
    end_SINR  = init_SINR+bin_size;
    
    bin_indexes = (EESM_vector>=init_SINR) & (EESM_vector<end_SINR);
    bin_SINRs   = EESM_vector(bin_indexes);
    bin_ACKs    = ACK_vector(bin_indexes);
    bin_sizes(bin_idx) = length(bin_SINRs);
    if bin_sizes(bin_idx)~=0
        bin_BLER(bin_idx)=1-mean(bin_ACKs);
    end
end