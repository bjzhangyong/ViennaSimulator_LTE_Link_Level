function [PDCCHmapping,UsedElements,NElements] = LTE_common_gen_PDCCH(Nrb,Nsc,Nsub,Ns,bandwidth,NoData)
% reserves space for the PDCCH (according to the test cases in TS 36.101 Appendix A3
% author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at


switch bandwidth
    case 1.4*10^6
        nsyms = 4;
    case {5*10^6,3*10^6}
        nsyms = 3;  
    case {10*10^6,15*10^6,20*10^6}
        nsyms = 2;   
end
PDCCHmapping = zeros(Nrb*Nsc,Nsub);
PDCCHmapping(:,1:nsyms) = 1;
PDCCHmapping(:,1:nsyms) = xor(PDCCHmapping(:,1:nsyms),NoData(:,1:nsyms));

UsedElements = zeros(Nrb*2,1);
for i2 = 1:2
    for i1 = 1:Nrb
        UsedElements(i1+(i2-1)*Nrb) = sum(sum(PDCCHmapping((i1-1)*Nsc+1:i1*Nsc,(i2-1)*Ns+1:Ns*i2)));
    end
end
NElements = sum(sum(UsedElements));   