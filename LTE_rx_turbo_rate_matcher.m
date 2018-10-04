function dk = LTE_rx_turbo_rate_matcher(LTE_params,ek,UE_signaling,UE,stream_index)
% LTE Turbo Code Rate Matcher, as of TS 36.212, Section 5.1.4.1.
% [dk] = LTE_rx_turbo_rate_matcher(ek,UE_signaling,UE)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   input_bits        ... ek (soft bits)
%           UE_signaling      ... BS signaling
% output:   output_bits       ... rate de-matched bits

% Nomenclature scheme
%        _______________________        ________________        ___________________________
%  dk-> | sub-block interleaver | vk-> | bit collection | wk-> | bit selection and pruning | ek->
%       |_______________________|      |________________       |___________________________|
%
% date of creation: 2008/08/11
% last changes: 
%   2008/09/15  Bosanska    added input [1 x 1]struct UE
%   2009/04/22  Jcolom      Changed how N_IR is set according to the new version of the standard (8.6.0)

C = length(ek);
wk = cell(1,C);
vk = cell(1,C);
dk = cell(1,C);
if UE.mode==3 || UE.mode==4
    K_MIMO = 2;
else
    K_MIMO = 1;
end
M_limit = 8;

N_IR = floor(UE.N_soft / (K_MIMO*min(LTE_params.HARQ_processes,M_limit)));

for i=1:C
    % Actual bit selection and pruning (genie driven)
    bit_selection_and_pruning_mapping = UE_signaling.turbo_rate_matcher(stream_index).bit_selection_and_pruning_mapping{i};
    wk{i} = LTE_common_turbo_rate_matching_bit_selection_and_pruning(ek{i},...
        bit_selection_and_pruning_mapping,...
        2,...
        UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).K_pi*3);
    
    % Circular buffer
    wk_length = length(wk{i});
    vk_length = wk_length/3;
    vk12 = reshape(wk{i}(vk_length+1:end),2,[]);
    vk{i} = [wk{i}(1:vk_length); vk12(1,:); vk12(2,:)];
    
    UE_signaling.turbo_rate_matcher(stream_index).ek_sizes(i) = length(bit_selection_and_pruning_mapping);
end
% Sub-block interleaver
dk = subblock_interleaver(LTE_params,UE_signaling.TB_segmentation(stream_index).C,UE_signaling.turbo_rate_matcher(stream_index),vk);
end

function dk = subblock_interleaver(LTE_params,C,signaling,vk)
% LTE Turbo Sub-Block deinterleaver, as of TS 36.212, Section 5.1.4.1.
% New implementation due to trying to find a deeply esoteric bug hidden
% somewhere in the channel coding and segmentation code.
% Author: Josep Colom Ikuno, josep.colom@nt.tuwien.ac.at
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2010/11/03
% last changes:

dk = cell(1,C);
for i_=1:C
    codeblock = vk{i_};
    R_tc  = signaling.subblock_interleaver(i_).R_tc;
    N_d   = signaling.subblock_interleaver(i_).Nd;
    d_k2  = LTE_common_soft_bit_interleaver(codeblock(3,:),signaling.subblock_interleaver(i_).vk2_mapping,0); % Interleave flag set to 0 (deinterleave)
    d_k01 = [
        reshape(codeblock(1,:),R_tc,[])
        reshape(codeblock(2,:),R_tc,[])
        ];
    d_k01 = d_k01(:,LTE_params.sub_block_interleaver_permutation_pattern_plus_one);
    d_k0  = reshape(d_k01(1:R_tc,:)',1,[]);
    d_k1  = reshape(d_k01((R_tc+1):end,:)',1,[]);
    % Final assembly removal of filler bits
    dk{i_} = [
        d_k0((N_d+1):end)
        d_k1((N_d+1):end)
        d_k2((N_d+1):end)
        ];
end

end
