function LLR_SD_C = LTE_softsphere(rx_layer_x,rx_user_symbols,Q,R,symbol_alphabet,bittable,nLayers,M)
% Soft Sphere decoder.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

siz = size(R,2);

if (siz < size(R,1)) % chop off unnecessary data
     R = R(1:siz,:);
     Q = Q(:,1:siz);
end

for ij = 1:size(rx_layer_x,2)
    for i = 1:nLayers
        [C,I] = min((abs(rx_layer_x(i,ij)*ones(1,2^M(i))-symbol_alphabet(i,1:2^M(i))).').^2);
        symbols_ZF(i,ij) = I.';    % ZF Symbols (integers)
    end
end
s_alph = [];
for mm=1:nLayers
    s_alph = [s_alph; symbol_alphabet(mm,symbols_ZF(mm,:))];
end

dist_ZF = sum(abs(Q'*rx_user_symbols.'-R*s_alph).^2,1);   % distance to the ZF solution initial value for SS Decoder
if (imag(R) == 0) % The SSD needs a complex matrix, or else the MEX version of it will crash
    R = complex(R);
end

% Soft Sphere Decoder
% Use this line if you would like to save the channel data for SSD testing
% save filename R Q rx_user_symbols dist_ZF symbols_ZF M symbol_alphabet bittable
LLR_SD_C = LTE_rx_soft_sd2(R,Q'*rx_user_symbols.',dist_ZF,int32(symbols_ZF),int32(M),symbol_alphabet.',bittable);
