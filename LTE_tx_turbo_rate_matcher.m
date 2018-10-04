function ek = LTE_tx_turbo_rate_matcher(LTE_params,dk,UE_signaling,UE,stream_index)
% LTE Turbo Code Rate Matcher, as of TS 36.212, Section 5.1.4.1.
% [ek UE_signaling] = LTE_tx_turbo_rate_matcher(dk,UE_signaling,UE)
% Author: Josep Colom Ikuno, josep.colom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   input_bits        ... Coded bits dk, as outputted by the turbo 
%                                 encoder This is a matrix of 3xN bits.
%           UE_signaling      ... BS signaling
% output:   output_bits       ... rate-matched bits
%           UE_signaling      ... UE_signaling with extra information, such as
%                                 number of padding bits used for
%                                 the interleaving (foe each sub-block).
%                                 Needed for the rate matching at the
%                                 receiver.

% Nomenclature scheme
%        _______________________        ________________        ___________________________
%  dk-> | sub-block interleaver | vk-> | bit collection | wk-> | bit selection and pruning | ek->
%       |_______________________|      |________________       |___________________________|
%
% date of creation: 2008/08/11
% last changes:
%   2009/04/22  Jcolom      Changed how N_IR is set according to the new version of the standard (8.6.0)

UE_signaling.turbo_rate_matcher(stream_index).total_bits_before_punturing = 0;
C = length(dk);
vk = cell(1,C);
wk = cell(1,C);
ek = cell(1,C);

% Sub-block interleaver
[   vk,... % cell array v_k
    UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver... % last entry (subblock_interleaver) is indexed with the CB index
    ] = subblock_interleaver(LTE_params,dk);

% Needed for the bit selection and transmission
if UE.mode==3 || UE.mode==4 || UE.mode==8 % According to 36.212, 5.1.4.1.2
    K_MIMO = 2;
else
    K_MIMO = 1;
end
M_limit = 8;
N_IR = floor(UE.N_soft / (K_MIMO*min(LTE_params.HARQ_processes,M_limit)));
K_pis = [UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver.K_pi];
K_w = 3*K_pis;
UE_signaling.turbo_rate_matcher(stream_index).total_bits_before_punturing = sum(K_w);

wk = cell(1,C);
ek = cell(1,C);

% Bit selection and pruning parameters
N_l    = UE_signaling.turbo_rate_matcher(stream_index).N_l;
Q_m    = UE_signaling.MCS_and_scheduling.CQI_params(stream_index).modulation_order;
G      = UE_signaling.turbo_rate_matcher(stream_index).G;
G_p    = G /(N_l*Q_m);
gamma  = mod(G_p,C);

% Necessary for balanced rate matching
TB_segmentation = UE_signaling.turbo_encoder(stream_index).encoded_CB_sizes; % UE_signaling.TB_segmentation(stream_index).CB_sizes;
ratios = TB_segmentation/sum(TB_segmentation);

if LTE_params.balanced_rate_matching && C > 1
    G_cb   = round(G*ratios);
    r_bits = G - sum(G_cb); % Maximum discrepancy of +- C bits
    if r_bits > 0
        M = mod(r_bits,C);
        N = floor(r_bits/C);
        CB_adjusting = N + [ones(1,M) zeros(1,C-M)];
    else
        M = mod(-r_bits,C);
        N = floor(-r_bits/C);
        CB_adjusting = -(N + [ones(1,M) zeros(1,C-M)]);
    end
    E_all = G_cb + CB_adjusting;
else
    G_p_C = G_p/C*ones(1,C);
    r = 0:(C-1);
    E_all = zeros(1,C);
    floored = (r<=(C-gamma-1));
    E_all(floored)  = N_l * Q_m * floor(G_p_C(floored));
    E_all(~floored) = N_l * Q_m * ceil(G_p_C(~floored));
end

for i=1:C
    UE_signaling.turbo_rate_matcher(stream_index).total_bits_before_punturing = UE_signaling.turbo_rate_matcher(stream_index).total_bits_before_punturing + 3*UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).K_pi;
    
    % Circular buffer
    wk{i} = [vk{i}(1,:) reshape(vk{i}(2:3,:),1,[])];
    
    % Bit selection and pruning parameters
    N_cb   = min(floor(N_IR/C),K_w(i));
    rv_idx = UE_signaling.turbo_rate_matcher(stream_index).rv_idx;
    K_pi   = K_pis(i);
    R_tc   = UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).R_tc;
    F      = UE_signaling.TB_segmentation(stream_index).F;
    r      = i-1;
    Nd     = UE_signaling.turbo_rate_matcher(stream_index).subblock_interleaver(i).Nd;
    
    E = E_all(i);
    k_0 = R_tc * (2*ceil(N_cb/(8*R_tc))*rv_idx+2);
    
    k_0_j = mod(k_0+(0:(2*E-1)),N_cb); % Assumed that there won't be so many filler bits
    if LTE_params.null_bit~=0
        null_pos       = wk{i}==LTE_params.null_bit; % The null positions in the codeblock
        not_null       = wk{i}(k_0_j+1)~=LTE_params.null_bit;
        k_0_j_not_null = k_0_j(not_null) + 1; % we store it in one-indexed form
    else
        null_pos       = [];
        k_0_j_not_null = k_0_j + 1; % we store it in one-indexed form
    end
    k_0_j_not_null = k_0_j_not_null(1:E);
    ek{i} = wk{i}(k_0_j_not_null);
    
    UE_signaling.turbo_rate_matcher(stream_index).k_0_j_not_null{i} = k_0_j_not_null; % the actual rate-matching mapping, such as  ek{i} = wk{i}(k_0_j_not_null);
    UE_signaling.turbo_rate_matcher(stream_index).ek_sizes(i)       = length(ek{i});
    UE_signaling.turbo_rate_matcher(stream_index).null_pos{i}       = null_pos; % The <NULL> positions in wk{i}
    UE_signaling.turbo_rate_matcher(stream_index).gamma             = gamma;
    UE_signaling.turbo_rate_matcher(stream_index).G_p               = G_p;
    
    % Backwards compatibility for the receiver
    UE_signaling.turbo_rate_matcher(stream_index).bit_selection_and_pruning_mapping{i} = k_0_j_not_null-1;
end
end

function [ v_k signaling ] = subblock_interleaver(LTE_params,dk)
% LTE Turbo Sub-Block interleaver, as of TS 36.212, Section 5.1.4.1.
% New implementation due to trying to find a deeply esoteric bug hidden
% somewhere in the channel coding and segmentation code.
% [y_perm Nd R K_pi] = LTE_common_subblock_interleaver(dk,interleave_flag,varargin)
% Author: Josep Colom Ikuno, josep.colom@nt.tuwien.ac.at
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2010/11/03
% last changes:

C = length(dk);
D = zeros(1,C);
for i_=1:C
    D(i_) = length(dk{i_});
end

C_tc = 32;
R_tc = ceil(D/C_tc);
N_d  = R_tc*C_tc - D;
K_pi = N_d + D;
Perm = LTE_params.sub_block_interleaver_permutation_pattern_plus_one;

% Fill in signaling information and perform colum permutation for first and
% second turbo-encoded rows.
v_k = cell(1,C);
R_start = 1;
for i_=1:C
    % d_k=0 and d_k=1
    this_CB = [LTE_params.null_bit*ones(3,N_d(i_),'uint8') dk{i_}];
    R_end   = R_start + 2*R_tc(i_) - 1;
    signaling(i_).Nd   = N_d(i_); % preserved this typo to ensure backwards compatibility
    signaling(i_).R_tc = R_tc(i_);
    signaling(i_).K_pi = K_pi(i_);
    v_k_01 = [
        reshape(this_CB(1,:)',C_tc,[])'
        reshape(this_CB(2,:)',C_tc,[])'
        ];
    v_k_01_perm = v_k_01(:,Perm);
    
    %d_k=2
    k = 0:K_pi(i_)-1;
    P = LTE_params.sub_block_interleaver_permutation_pattern(floor(k/R_tc(i_))+1); % From TS 36.212, Table 5.1.4-1
    k_mod_R = mod(k,R_tc(i_));
    interleaving_map = mod(P+C_tc*k_mod_R+1,K_pi(i_));
    v_k2 = LTE_common_bit_interleaver(this_CB(3,:),interleaving_map,1); % interleave_flag=1
    signaling(i_).vk2_mapping = interleaving_map;
    
    % Assemble permutated codeblock
    v_k{i_} = [
        reshape(v_k_01_perm(1:R_tc(i_),:),1,[])
        reshape(v_k_01_perm(R_tc(i_)+1:end,:),1,[])
        v_k2
        ];
end

end
