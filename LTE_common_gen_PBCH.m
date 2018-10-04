function [PBCHmapping,UsedElements,NElements,NoData_indices] = LTE_common_gen_PBCH(Nrb,Nsc,Nsub,NIDcell,subframe_num,Ns,nTX)
% reserves space for the PBCH
% author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at

[RefSym,RefMapping,total_no_refsym,NoData_indices] = LTE_common_gen_Reference_Signal(Nrb, Nsub, nTX, NIDcell, subframe_num); 

if subframe_num ~= 1
    PBCHmapping = zeros(size(NoData_indices));
    NElements = 0;
    UsedElements = zeros(Nrb*2,1);
else
%     if strcmp(CyclicPrefix,'normal')
%         nbits = 1920;
%     else
%         nbits = 1728;
%     end
    
    k = Nrb*Nsc/2-36;
    PBCHmapping = zeros(size(NoData_indices));
    for i1 = 0:71
        for i2 = 0:3
            if ~NoData_indices(k+i1+1,i2+8)
               PBCHmapping(k+i1+1,i2+8) = 1;
            end
        end
    end
       
    UsedElements = zeros(Nrb*2,1);
    for i2 = 1:2
        for i1 = 1:Nrb
            UsedElements(i1+(i2-1)*Nrb) = sum(sum(PBCHmapping((i1-1)*Nsc+1:i1*Nsc,(i2-1)*Ns+1:Ns*i2)));
        end
    end
    NElements = sum(sum(UsedElements));    
end

