function [RefSym, RefMapping, no_refsym_per_slot, NoData_indices] = LTE_common_gen_Reference_Signal(Nrb, Nsub, nTX, NIDcell, subframe)
% Reference signal generation and mapping.
% Generation and Mapping based on one time subframe of 1 ms - 3GPP TS 36.211 V.8.2.0 section 6.10,
% [RefSym, RefMapping] = LTE_Common_gen_Reference_Signal(Nrb, Ns, AtPort, NIDcell, subframe)
% Author: Qi Wang, qwang@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   Nrb                 ... [1 x 1]double number of resource blocks
%           Nsub                ... [1 x 1]double number of OFDM symbols in one subframe
%           nTX                 ... [1 x 1]double number of tranmit antenna
%           NIDcell             ... [1 x 1]double cell identity
%           subframe            ... [1 x 1]double number of the subframe transmitted             
% output:   RefSym              ... [m x n]double generated reference signals
%           RefMapping          ... [LTE_params.Nsc*Nrb x Nsub x nTX]logical mapping matrix
%           no_refsym_per_slot  ... [1 x 1]double number of resource elements occupied by reference symbols per slot
%           NoData_indices      ... [LTE_params.Nsc*Nrb x Nsub]logical matrix for non-data symbols
%
% date of creation: 2008/08/11
% last changes:
%             2009/03/02    qwang   adjust to MIMO

AtPort_vec = 0:nTX-1;
RefSym = zeros(2*Nrb,4,nTX);
RefMapping = false(12*Nrb, Nsub, nTX);
NoData_indices = 0;
NmaxRB = 110;

% Number of resource elements occupied by reference symbols per slot
switch nTX
    case 1
        no_refsym_per_slot = 4;
    case 2
        no_refsym_per_slot = 8;
    case 4
        no_refsym_per_slot = 12;
    otherwise
        error('wrong nTX.')
end

% Generation
NumBit = NmaxRB*4;
% Pseudo-random sequence sec:7.2
pn_gen_x1 = commsrc.pn('GenPoly', [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1], ...
              'InitialStates',   [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1],   ...
              'Shift',         1600,               ...
              'NumBitsOut',    NumBit);
pn_seq1 = generate(pn_gen_x1);

for aa = 1:nTX
    AtPort = AtPort_vec(aa);
    clear RefSym_temp RefMapping_temp;

    % initial sequence for x2, TS36.211 V8.2.0 6.10.1.1
    switch AtPort
        case {0,1}
            SyminSub_i = [1, Nsub/2-2, Nsub/2+1, Nsub-2];
        case {2,3}
            SyminSub_i = [2, Nsub/2+2];
        otherwise 
            error('Reference signal for this antenna port is not implemented.');
    end
    c_ini = de2bi(2^10*(7*(subframe*2-1+floor(SyminSub_i/Nsub*2))+mod(SyminSub_i,Nsub/2)+1)*(2*NIDcell+1)+2*NIDcell+(Nsub==14),31);
    % old version 36.211 8.2.0
%     c_ini = de2bi(2^13*SyminSub_i+2^9*(subframe-1)+NIDcell,31);
    pn_seq2 = nan(NumBit,length(SyminSub_i));
    pn_seq = nan(NumBit,length(SyminSub_i));

    for ii=1:size(c_ini,1)
        pn_gen_x2 = commsrc.pn('GenPoly', [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1], ...
                      'InitialStates',   c_ini(ii,:),   ...
                      'Shift',         1600,               ...
                      'NumBitsOut',    NumBit);
        pn_seq2(:,ii) = generate(pn_gen_x2);
        pn_seq(:,ii) = mod((pn_seq1+pn_seq2(:,ii)),2);
    end

    r = 1/sqrt(2)*(1-2*pn_seq(1:2:end,:))+1i*1/sqrt(2)*(1-2*pn_seq(2:2:end,:));
    RefSym_temp = r((1:Nrb*2)+NmaxRB-Nrb,:);
    RefSym(1:size(RefSym_temp,1),1:size(RefSym_temp,2),aa) = RefSym_temp;

    % Mapping
    RefMapping_temp = false(12*Nrb, Nsub);

    v_shift = mod(NIDcell,6);
    k = zeros(Nrb*2,length(SyminSub_i));
    switch AtPort
        case 0
            k(:,[1 3]) = repmat(6*(0:Nrb*2-1)+1 + v_shift,[2 1])';
            k(:,[2 4]) = repmat(6*(0:Nrb*2-1)+1 + mod((3+v_shift),6),[2 1])';
        case 1
            k(:,[1 3]) = repmat(6*(0:Nrb*2-1)+1 + mod((3+v_shift),6),[2 1])';
            k(:,[2 4]) = repmat(6*(0:Nrb*2-1)+1 + v_shift,[2 1])';
        case 2
            k(:,1) = 6*(0:Nrb*2-1)+1 + v_shift;
            k(:,2) = 6*(0:Nrb*2-1)+1 + mod((3+v_shift),6);
        case 3
            k(:,1) = 6*(0:Nrb*2-1)+1 + mod((3+v_shift),6);
            k(:,2) = 6*(0:Nrb*2-1)+1 + v_shift;
        otherwise
            error('Mapping for this Antenna Port is not implemented.');
    end
    l = SyminSub_i;

    for ii=1:length(l)
        RefMapping_temp(k(:,ii),l(ii))=1;
    end
    RefMapping(:,:,aa) = RefMapping_temp;
    
    NoData_indices = RefMapping(:,:,aa)|NoData_indices;
end
