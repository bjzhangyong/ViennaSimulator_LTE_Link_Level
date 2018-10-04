function [precode_y PRE] = LTE_precoding(tx_mode, layer_x, nAtPort, codebook_index, nLayers,indices, Ntot,CDD,LTE_params,slots,V,U)
% This function does the LTE Precoding for the different transmission modes
% according to standards 3GPPP TS 36.211 and TS 36.213
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at
switch tx_mode
    case 1 % SISO
        precode_y = layer_x;
    %     W = 1;          % just necessary to have a receiver which is the same for SISO and MIMO
    %     U = 1;
    %     D = 1;
        PRE = 1;
       
    case 2 % Transmit Diversity
        [Z W U D] = LTE_common_get_precoding_matrix(tx_mode,nAtPort,codebook_index,nLayers,LTE_params);
        X = 1/sqrt(2)*Z*[real(layer_x(:,:));imag(layer_x(:,:))];
        c = length(layer_x(1,:));
        precode_y = zeros(nAtPort,nAtPort*c);
        for i = 1:nAtPort
            precode_y(:,i:nAtPort:nAtPort*c-(nAtPort-i))=X((i-1)*nAtPort+1:i*nAtPort,:);
        end
        PRE = 1;

    case 3 % Open loop Spatial Multiplexing (uses large delay CDD according to TS 36.213 section 7.1.3)
       [Z W U D] = LTE_common_get_precoding_matrix(tx_mode,nAtPort,codebook_index,nLayers,LTE_params);  
       if (nAtPort == 2)    % faster encoding (for 2 Antennaports W is not dependent on the time index)
           PRE = zeros(nAtPort,nLayers,2);
           PRE(:,:,1) = W*D*U;
           PRE(:,:,2) = W*D^2*U;
           precodey1= PRE(:,:,1)*layer_x(:,1:2:end);
           precodey2= PRE(:,:,2)*layer_x(:,2:2:end);
           precode_y = reshape([precodey1;precodey2],nAtPort,[]);
       else
       l = 1:length(layer_x);
       k = mod(floor((l-1)/nLayers),4)+1;
       p = mod(l-1,nLayers)+1;
       period = lcm(seqperiod(k),seqperiod(p));
       PRE = zeros(nAtPort,nLayers,period);
       precode_y = zeros(nAtPort,length(layer_x));
       for i = 1:period
           PRE(:,:,i)= W(:,:,k(i))*D^(p(i))*U;
           temp = PRE(:,:,i)*layer_x(:,i:period:end);
           precode_y(:,i:period:size(temp,2)*period) = temp;
       end
       end
   
    case 4 % Closed loop Spatial Multiplexing
    [ind1,ind2] = size(codebook_index);
    W = zeros(nAtPort,nLayers,ind1,ind2);
    for i1 = 1:ind1
        for i2 = 1:ind2
            [Z W(:,:,i1,i2) U D] = LTE_common_get_precoding_matrix(tx_mode,nAtPort,codebook_index(i1,i2),nLayers,LTE_params);          
        end
    end
    precode_y = zeros(nAtPort,length(layer_x));
    PRE = zeros(nAtPort,nLayers,2,Ntot); % the 2 is for the two slots per subframe
%     if (CDD(1) == 0)
    for i1 = 1:ind1
        for i2 = 1:ind2
            PRE(:,:,i2,(LTE_params.Nsc)*(i1-1)+1:LTE_params.Nsc*i1) = repmat(W(:,:,i1,i2),[1,1,1,LTE_params.Nsc]);
        end
    end

%         precode_y = PRE(:,:,1)*layer_x;
%     else
%     for i = 1:length(layer_x);
%            precode_y(:,i) = D^(indices(i))*W*layer_x(:,i);
%     end
     for i1 = 1:Ntot
         for i2 = 1:2
            ind1 = ~(indices-i1);
            ind2 = ~(slots-i2);
            ind = ind1 & ind2;
            if CDD(1) ~= 0
                PRE(:,:,i2,i1)=D^(i1)*PRE(:,:,i2,i1);
            end
            precode_y(:,ind) =PRE(:,:,i2,i1)*layer_x(:,ind);
         end
     end
%     end

    case 6  % Interference Alignment
        for i = 1:size(layer_x,2)
            if LTE_params.IA_time_granularity == 14
                time_index = 1;
            elseif LTE_params.IA_time_granularity == 7
                time_index = LTE_params.IA_time_granularity*(slots(i)-1)+1;
            else
                time_index = 1; % at the moment only slot (7) and subframe (14) time granularity is implemented
            end

            precode_y(:,i) = V{indices(i),time_index}*layer_x(:,i);
        end
        PRE = eye(4,2);
        
    otherwise
        error('TX mode %d not supported.',tx_mode);
end
