function [LLR_SD,M] = LTE_detect_SISO(MCS_and_scheduling,filtering,rx_user_symbols,H_est,H_est_user,LTE_params,receiver)
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

if (strcmp(filtering,'BlockFading'))
     for ctr = 1:LTE_params.Ntot
            ind = ~(indices-ctr);
            if (~ind)
                continue
            end
            H_complete = squeeze(H_est(ctr,1,:));
            rx_layer_x = 1/H_complete*rx_user_symbols(ind,:).';
            switch receiver
                case 'SSD'
                    [Q,R] = qr(H_complete);
                    LLR_SD(:,ind) = LTE_softsphere(rx_layer_x,rx_user_symbols(ind,:),Q,R,symbol_alphabet,bittable,nLayers,M);
                case 'ZF'
                     LLR_SD(:,ind) = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M);
            end
     end
else
    for ctr = 1:l(end)
            H_complete = squeeze(H_est_user(ctr,:));
            rx_layer_x = 1/H_complete*rx_user_symbols(ctr,:).';
            switch receiver
                case 'SSD'
                    [Q,R] = qr(H_complete);
                    LLR_SD(:,ctr) = LTE_softsphere(rx_layer_x,rx_user_symbols(ctr,:),Q,R,symbol_alphabet,bittable,nLayers,M);
                case 'ZF'
                     LLR_SD(:,ctr) = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M);
            end
    end
end