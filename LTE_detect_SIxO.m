function [LLR_SD,M,rx_layer_x_equalized] = LTE_detect_SIxO(MCS_and_scheduling,filtering,rx_user_symbols,H_est,H_est_user,LTE_params,receiver,receiver_k,sigma_n2,MSE)
% SISO detection.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

nLayers = MCS_and_scheduling.nLayers;
M = MCS_and_scheduling.CQI_params.modulation_order;
bittable = false(M,2^max(M));
symbol_alphabet = zeros(1,2^max(M));
bittable(1:M,1:2^M)=LTE_params.bittable{M}; % Bitmapping table
symbol_alphabet(1,1:2^M)=LTE_params.SymbolAlphabet{M}.'; % Symbol alphabet

LLR_SD = zeros(M,length(rx_user_symbols));   % Log likelihood Ratios of the Spere decoder

l = 1:length(rx_user_symbols);
indices = MCS_and_scheduling.freq_indices;

rx_layer_x_equalized = zeros(length(rx_user_symbols),nLayers);

if (strcmp(filtering,'BlockFading'))
     for ctr = 1:LTE_params.Ntot
            ind = ~(indices-ctr);
            if (~ind)
                continue
            end
            H_complete = squeeze(H_est(ctr,1,:));
            switch receiver
                case 'SSD'
                    rx_layer_x = pinv(H_complete)*rx_user_symbols(ind,:).';
                    [Q,R] = qr(H_complete);
                    LLR_SD(:,ind) = LTE_softsphere(rx_layer_x,rx_user_symbols(ind,:),Q,R,symbol_alphabet,bittable,nLayers,M);
                case 'SSDKB'
                    rx_layer_x = pinv(H_complete)*rx_user_symbols(ind,:).';
                    [Q,R] = qr(H_complete);
                    LLR_SD(:,ind) = LTE_softsphere_kbest(rx_layer_x,rx_user_symbols(ind,:),Q,R,symbol_alphabet,bittable,nLayers,M,receiver_k);
				case 'ZF'
                    inv_temp = pinv(H_complete);
                    rx_layer_x = inv_temp*rx_user_symbols(ind,:).';
                    rx_layer_x_equalized(ind,:) = rx_layer_x;
                    Hg = inv_temp*H_complete;
                    noise_enhancement = sum(abs(inv_temp).^2,2)*ones(M,length(rx_layer_x));
                    LLR_SD(:,ind) = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                case 'MMSE'
                    temp = H_complete'*H_complete;
%                     inv_temp = (temp+(sigma_n2+MSE(ctr))*eye(size(temp)))^-1*H_complete'; 
                    inv_temp = (temp+sigma_n2*eye(size(temp)))^-1*H_complete';
                    rx_layer_x = inv_temp*rx_user_symbols(ind,:).';
                    rx_layer_x_equalized(ind,:) = rx_layer_x;
                    Hg = inv_temp*H_complete;
                    noise_enhancement = sum(abs(inv_temp).^2,2)*ones(M,length(rx_layer_x));
                    LLR_SD(:,ind) = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
            end
     end
else
    for ctr = 1:l(end)
            H_complete = squeeze(H_est_user(ctr,:)).';
            switch receiver
                case 'SSD'
                    rx_layer_x = pinv(H_complete)*rx_user_symbols(ctr,:).';
                    [Q,R] = qr(H_complete);
                    LLR_SD(:,ctr) = LTE_softsphere(rx_layer_x,rx_user_symbols(ctr,:),Q,R,symbol_alphabet,bittable,nLayers,M);
                case 'SSDKB'
                    rx_layer_x = pinv(H_complete)*rx_user_symbols(ctr,:).';
                    [Q,R] = qr(H_complete);
                    LLR_SD(:,ctr) = LTE_softsphere_kbest(rx_layer_x,rx_user_symbols(ctr,:),Q,R,symbol_alphabet,bittable,nLayers,M,receiver_k);
				case 'ZF'
                    inv_temp = pinv(H_complete);
                    rx_layer_x = inv_temp*rx_user_symbols(ctr,:).';
                    Hg = inv_temp*H_complete;
                    noise_enhancement = sum(abs(inv_temp).^2,2)*ones(M,length(rx_layer_x));
                    LLR_SD(:,ctr) = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                case 'MMSE'
                    temp = H_complete'*H_complete;
                    inv_temp = (temp+sigma_n2*eye(size(temp)))^-1*H_complete';
                    rx_layer_x = inv_temp*rx_user_symbols(ctr,:).';
                    Hg = inv_temp*H_complete;
                    noise_enhancement = sum(abs(inv_temp).^2,2)*ones(M,length(rx_layer_x));
                    LLR_SD(:,ctr) = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);                    
            end
    end
end
