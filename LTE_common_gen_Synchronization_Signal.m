function[PrimSync, PrimMapping, SecSync, SecMapping, SyncUsedElements,ReservedMapping,ResUsedElements,ResNElements] = LTE_common_gen_Synchronization_Signal(NIDcell, Nrb, Nsc,Ns,subframe)
% Synchronization signal generation and mapping
% Generation and Mapping based on one subframe slot of 1 ms, 3GPP TS 36.211 V.8.2.0 section 6.11
% [PrimSync, PrimMapping, SecSync, SecMapping] = LTE_Common_gen_Synchronization_Signal(NIDcell, Nrb, Nsc,Ns,subframe)
% Author: Qi Wang, qwang@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   NIDcell             ... [1 x 1]double cell identity
%           Nrb                 ... [1 x 1]double number of resource blocks
%           Nsc                 ... [1 x 1]double number of subcarriers in one resource block
%           Ns                  ... [1 x 1]double number of OFDM symbols in one slot
%           subframe            ... [1 x 1]double number of the subframe transmitted             
% output:   PrimSync            ... [62 x 1]double generated primary synchronization signals
%           PrimMapping         ... [LTE_params.Nsc*Nrb x 2*Ns]logical mapping matrix
%           SecSync             ... [62 x 1]double generated primary synchronization signals
%           SecMapping          ... [LTE_params.Nsc*Nrb x 2*Ns]logical mapping matrix
%
% date of creation: 2008/08/11
% last changes: 2008/09/03  Bosanska    added output SyncUsedElements for multi-
%                                       -user implementation and resource mapping

NID1 = floor(NIDcell/3);
NID2 = mod(NIDcell,3);

%% Primary synchronization signal
switch NID2
    case 0 
        u = 25;
    case 1
        u = 29;
    case 2
        u = 34;
    otherwise
         error('invalid NID2.');
end
PrimSync_temp = zeros(31,2);
for n=1:31
    PrimSync_temp(n,1) = exp(-1i*pi*u*(n-1)*n/63);
    PrimSync_temp(n,2) = exp(-1i*pi*u*(n+31)*(n+32)/63);
end

PrimSync = reshape(PrimSync_temp,62,[]);

% Primary synchronization signal mapping and used resource elements per
% resource blocks
PrimMapping = false(Nrb*Nsc,Ns);
SyncUsedElements = zeros(Nrb,2);
k = zeros(62,1);
for n=1:62
    k(n)=n-32+Nrb*Nsc/2;
end
k = k+1;    % corrected to matlab indices (zero excluded) - without this the Synchronisationsignal is not exactly in the middle
PrimMapping(k,Ns/2)=1;

SyncUsedElements(ceil(k(1)/Nsc),1) = 2*(ceil(k(1)/Nsc)*Nsc-k(1)+1);
% SyncUsedElements(ceil(k(1)/Nsc),1) = 2*(ceil(k(1)/Nsc)*Nsc-k(1)+1);
if(mod(k(end),Nsc) == 0)
    SyncUsedElements(ceil(k(1)/Nsc)+1:ceil(k(end)/Nsc),1) = 2*12;
else
    SyncUsedElements(ceil(k(1)/Nsc)+1:ceil(k(end)/Nsc)-1,1) = 2*12;
    SyncUsedElements(ceil(k(end)/Nsc),1) = 2*12-2*(ceil(k(end)/Nsc)*Nsc-k(end));
%     SyncUsedElements(ceil(k(end)/Nsc),1) = 2*(ceil(k(end)/Nsc)*Nsc-k(end));
end

%% Secondary synchronization signal
qq = floor(NID1/30);
q = floor((NID1+qq*(qq+1)/2)/30);
mm = NID1+q*(q+1)/2;
m0 = mod(mm,31);
m1 = mod(m0+floor(mm/31)+1, 31);

PN_gen = commsrc.pn('GenPoly',       [1 0 0 1 0 1], ...
                     'InitialStates', [1 0 0 0 0],   ...
                     'Shift',         0,   ...
                     'NumBitsOut',    31);
PN_seq_temp = generate(PN_gen).';
PN_seq = 1-2*PN_seq_temp;

SC_gen1 = commsrc.pn('GenPoly',       [1 0 1 0 0 1], ...
                     'InitialStates', [1 0 0 0 0],   ...
                     'Shift',         0,   ...
                     'NumBitsOut',    31);
SC_seq1_temp = generate(SC_gen1).';
SC_seq1 = 1-2*SC_seq1_temp;

SC_gen2 = commsrc.pn('GenPoly',       [1 1 0 1 1 1], ...
                     'InitialStates', [1 0 0 0 0],   ...
                     'Shift',         0,   ...
                     'NumBitsOut',    31);
SC_seq2_temp = generate(SC_gen2).';
SC_seq2 = 1-2*SC_seq2_temp;

SecSync = zeros(62,1);

for n=0:30
    switch subframe
        case 1
            SecSync(2*n+1) = PN_seq(mod(n+m0,31)+1).*SC_seq1(mod(n+NID2, 31)+1);
            SecSync(2*n+2) = PN_seq(mod(n+m1,31)+1).*SC_seq1(mod(n+NID2+3,31)+1).*SC_seq2(mod(n+mod(m0,8),31)+1);
        case 6
            SecSync(2*n+1) = PN_seq(mod(n+m1,31)+1).*SC_seq1(mod(n+NID2, 31)+1);
            SecSync(2*n+2) = PN_seq(mod(n+m0,31)+1).*SC_seq1(mod(n+NID2+3,31)+1).*SC_seq2(mod(n+mod(m1,8),31)+1);
        otherwise     
    end
end

% Secondary synchronization signal mapping
SecMapping = false(Nrb*Nsc,Ns);
k = zeros(62,1);
for n=1:62
    k(n)=n-32+Nrb*Nsc/2;
end
k = k+1;
SecMapping(k,Ns/2-1)=1;
ReservedMapping = zeros(size(SecMapping));
ReservedMapping([k(1)-5:k(1)-1,k(end)+1:k(end)+5],Ns/2-1:Ns/2) = 1;  % this space is reserved to keep ICI low
ResUsedElements = zeros(Nrb*2,1);
for i2 = 1:2
    for i1 = 1:Nrb
        ResUsedElements(i1+(i2-1)*Nrb) = sum(sum(ReservedMapping((i1-1)*Nsc+1:i1*Nsc,(i2-1)*Ns/2+1:Ns/2*i2)));
    end
end
ResNElements = sum(sum(ResUsedElements));   
