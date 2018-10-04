function [Z W U D] = LTE_common_get_precoding_matrix(tx_mode,nAtPort,codebook_index,nLayers,LTE_params)
% This function returns the precoding matrix for a specific transmission mode
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

Z = 0;
W = 0;
D = 0;
U = 0;  % the matrices which are not used for the specific tx_mode are set to zero

if (tx_mode == 2)   % transmit diversity
    Z = LTE_params.Z{nAtPort/2};    % matrix corresponding to 36.211 section 6.3.4.3  
elseif (tx_mode == 3)   % open loop spatial multiplexing section 6.3.4.2.2
    U = LTE_params.U_l{nLayers};
    D = LTE_params.D_l{nLayers};
    if (nAtPort == 2)
        W = zeros(2,nLayers,numel(codebook_index));
        W(:,:) = LTE_params.W{nLayers}(:,:,codebook_index+1);
    else
        W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
        W = zeros(4,nLayers,4);
        for ii = 13:16
            W(:,:,ii-12) = W_temp(:,LTE_params.mapping{nLayers}(ii,:),ii-12);
        end
    end
elseif (tx_mode == 4)    % closed loop spatial multiplexing section 6.3.4.2.1
    D = LTE_params.D_s{nAtPort/2};
    if (nAtPort == 2)
        W = LTE_params.W{nLayers}(:,:,codebook_index+1);
    else
        W_temp = 1/sqrt(nLayers)*LTE_params.Wn(:,:,codebook_index+1);
        W = W_temp(:,LTE_params.mapping{nLayers}(codebook_index+1,:),1);
    end
end