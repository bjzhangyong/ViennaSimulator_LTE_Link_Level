function LLR_SD_C = LTE_demapper(rx_layer_x,symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement)
% Soft Sphere decoder.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

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
dist_ZF = sum(abs(rx_layer_x-s_alph).^2,1);   % distance to the ZF solution initial value for SS Decoder

% Soft Sphere Decoder
if imag(Hg) == 0
    Hg = complex(Hg);
end

LLR_SD_C = LTE_rx_soft_sd2(Hg,rx_layer_x,dist_ZF,int32(symbols_ZF),int32(M),symbol_alphabet.',bittable)./noise_enhancement;



%% If you use ZF and there is no interference (inter carrier, imperfect channel knowledge,...) use this demapper to increase the speed (it demaps the layers independently)
% for ij = 1:size(rx_layer_x,2)
%     for i = 1:nLayers
%         [C,I] = min((abs(rx_layer_x(i,ij)*ones(1,2^M(i))-symbol_alphabet(i,1:2^M(i))).').^2);
%         symbols_ZF(i,ij) = I.';    % ZF Symbols (integers)
%     end
% end
% s_alph = [];
% for mm=1:nLayers
%     s_alph = [s_alph; symbol_alphabet(mm,symbols_ZF(mm,:))];
% end
% % dist_ZF = sum(abs(rx_layer_x-s_alph).^2,1);   % distance to the ZF solution initial value for SS Decoder
% dist_ZF = abs(rx_layer_x-s_alph).^2;   % distance to the ZF solution initial value for SS Decoder
% % Soft Sphere Decoder
% % if imag(Hg) == 0
% %     Hg = complex(Hg);
% % end
% 
% LLR_SD_C = zeros(size(noise_enhancement));
% summe = 1;
% for i = 1:nLayers
%     LLR_SD_C(summe:summe+M(i)-1,:) = LTE_rx_soft_sd2(1+eps*1i,rx_layer_x(i,:),dist_ZF(i,:),int32(symbols_ZF(i,:)),int32(M(i)),symbol_alphabet(i,:).',bittable(summe:summe+M(i)-1,:))./noise_enhancement(summe:summe+M(i)-1,:);
%     summe = summe+M(i);
% end

