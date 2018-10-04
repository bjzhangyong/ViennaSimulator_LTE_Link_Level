function [LLR_SD,M,rx_layer_x_equalized] = LTE_detect_TxD(MCS_and_scheduling,BS_nAtPort,rx_user_symbols,UE_nRX,H_est_user,LTE_params,receiver,receiver_k,sigma_n2)
% Transmit diversity detection.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at



if (BS_nAtPort == 2)
    nLayers = 2;
    rx_layer_x_equalized = zeros(length(rx_user_symbols)/2,nLayers);
    
    M = [MCS_and_scheduling.CQI_params.modulation_order] * ones(1,nLayers);    % Modulation order
    bittable = false(sum(M(1:nLayers)),2^max(M));
    symbol_alphabet = zeros(nLayers,2^max(M));
    for i = 1:nLayers
        bittable(sum(M(1:i-1))+(1:M(i)),1:2^M(i))=LTE_params.bittable{M(i)}; % Bitmapping table
        symbol_alphabet(i,1:2^M(i))=LTE_params.SymbolAlphabet{M(i)}.'; % Symbol alphabet
    end
    symbol_alphabet(2,:) = conj(symbol_alphabet(2,:));  % this is because also the received symbols get conjugated
    rx_layer_x = zeros(2,length(rx_user_symbols)/2);    % initialization of variable for received layers
    LLR_SD = zeros(M(1)*nLayers,length(rx_user_symbols)/2);   % Log likelihood Ratios of the Spere decoder
    
    if (UE_nRX == 2)                        % 2x2 case
        for i = 1:2:length(rx_user_symbols)-1
            rx_user_symbols_combined = [rx_user_symbols(i,1);rx_user_symbols(i,2);conj(rx_user_symbols(i+1,1));conj(rx_user_symbols(i+1,2))];
            H_est_user_combined = [H_est_user(i,1,1),-H_est_user(i,1,2);H_est_user(i,2,1),-H_est_user(i,2,2);conj(H_est_user(i+1,1,2)),conj(H_est_user(i+1,1,1));conj(H_est_user(i+1,2,2)),conj(H_est_user(i+1,2,1))];  % generate equivalent channel matrix
            switch receiver
                case 'SSD'
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = pinv(H_est_user_combined)*rx_user_symbols_combined*sqrt(2);   % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                      % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:);
                    %%%%%%%%%%%%%%%%%%%
                    [Q,R]= qr(H_est_user_combined);
                    LLR_SD(:,(i+1)/2) = LTE_softsphere(rx_layer_x(:,(i+1)/2),rx_user_symbols_combined.'*sqrt(2),Q,R,symbol_alphabet,bittable,nLayers,M);
                case 'SSDKB'
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = pinv(H_est_user_combined)*rx_user_symbols_combined*sqrt(2);   % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                      % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:);
                    %%%%%%%%%%%%%%%%%%%
                    [Q,R]= qr(H_est_user_combined);
                    LLR_SD(:,(i+1)/2) = LTE_softsphere_kbest(rx_layer_x(:,(i+1)/2),rx_user_symbols_combined.'*sqrt(2),Q,R,symbol_alphabet,bittable,nLayers,M,receiver_k);
				case 'ZF'
                    inv_temp = pinv(H_est_user_combined);
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);   % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                      % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x_equalized((i+1)/2,1) = rx_user_symbols_combined_2(1,:); 
                    rx_layer_x_equalized((i+1)/2,2) = rx_user_symbols_combined_2(2,:)';
                    %%%%%%%%%%%%%%%%%%%
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1)];
                    LLR_SD(:,(i+1)/2) = LTE_demapper(rx_layer_x(:,(i+1)/2),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                case 'MMSE'
                    temp = H_est_user_combined'*H_est_user_combined;
                    inv_temp = (temp+sigma_n2*eye(size(temp)))^-1*H_est_user_combined';
                    % MMSE Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+1)/2,1) = rx_user_symbols_combined_2(1,:); 
                    rx_layer_x_equalized((i+1)/2,2) = rx_user_symbols_combined_2(2,:);
                    %%%%%%%%%%%%%%%%%%%
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1)];
                    LLR_SD(:,(i+1)/2) = LTE_demapper(rx_layer_x(:,(i+1)/2),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
            end
            
        end
        
    elseif (UE_nRX == 1)                % 2x1 case
        for i = 1:2:length(rx_user_symbols)-1
            rx_user_symbols_combined = [rx_user_symbols(i,1);conj(rx_user_symbols(i+1,1))];
            H_est_user_combined = [H_est_user(i,1,1),-H_est_user(i,1,2);conj(H_est_user(i+1,1,2)),conj(H_est_user(i+1,1,1))]; % generate equivalent channel matrix
            
            switch receiver
                case 'SSD'
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = pinv(H_est_user_combined)*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:);
                    %%%%%%%%%%%%%%%%%%%
                    [Q,R]= qr(H_est_user_combined);
                    LLR_SD(:,(i+1)/2) = LTE_softsphere(rx_layer_x(:,(i+1)/2),rx_user_symbols_combined.'*sqrt(2),Q,R,symbol_alphabet,bittable,nLayers,M);
                case 'SSDKB'
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = pinv(H_est_user_combined)*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:);
                    %%%%%%%%%%%%%%%%%%%
                    [Q,R]= qr(H_est_user_combined);
                    LLR_SD(:,(i+1)/2) = LTE_softsphere_kbest(rx_layer_x(:,(i+1)/2),rx_user_symbols_combined.'*sqrt(2),Q,R,symbol_alphabet,bittable,nLayers,M,receiver_k);
				case 'ZF'
                    inv_temp = pinv(H_est_user_combined);
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x_equalized((i+1)/2,1) = rx_user_symbols_combined_2(1,:); 
                    rx_layer_x_equalized((i+1)/2,2) = rx_user_symbols_combined_2(2,:)';
                    %%%%%%%%%%%%%%%%%%%
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1)];
                    LLR_SD(:,(i+1)/2) = LTE_demapper(rx_layer_x(:,(i+1)/2),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                case 'MMSE'
                    temp = H_est_user_combined'*H_est_user_combined;
                    inv_temp = (temp+sigma_n2*eye(size(temp)))^-1*H_est_user_combined';
                    % MMSE Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+1)/2)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+1)/2)=rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+1)/2,1) = rx_user_symbols_combined_2(1,:); 
                    rx_layer_x_equalized((i+1)/2,2) = rx_user_symbols_combined_2(2,:);
                    %%%%%%%%%%%%%%%%%%%
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1)];
                    LLR_SD(:,(i+1)/2) = LTE_demapper(rx_layer_x(:,(i+1)/2),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                    
            end
        end
    else error('number of RX antennas not supported');
    end
else
    nLayers = 4;
    M = [MCS_and_scheduling.CQI_params.modulation_order] * ones(1,nLayers);    % Modulation order
    rx_layer_x_equalized = zeros(length(rx_user_symbols)/4,nLayers);
    bittable = false(sum(M(1:nLayers)),2^max(M));
    symbol_alphabet = zeros(nLayers,2^max(M));
    for i = 1:nLayers
        bittable(sum(M(1:i-1))+(1:M(i)),1:2^M(i))=LTE_params.bittable{M(i)}; % Bitmapping table
        symbol_alphabet(i,1:2^M(i))=LTE_params.SymbolAlphabet{M(i)}.'; % Symbol alphabet
    end
    symbol_alphabet(2,:) = conj(symbol_alphabet(2,:));  % this is because also the received symbols get conjugated
    symbol_alphabet(4,:) = conj(symbol_alphabet(4,:));  % this is because also the received symbols get conjugated
    LLR_SD = zeros(M(1)*nLayers,length(rx_user_symbols)/4);   % Log likelihood Ratios of the Spere decoder
    rx_layer_x = zeros(4,length(rx_user_symbols)/4);    % initialization of variable for received layers

    if (UE_nRX == 4)                % 4x4 case
        for i = 1:4:length(rx_user_symbols)-3
            rx_user_symbols_combined = [rx_user_symbols(i,1);rx_user_symbols(i,2);conj(rx_user_symbols(i+1,1));conj(rx_user_symbols(i+1,2));
                rx_user_symbols(i+2,1);rx_user_symbols(i+2,2);conj(rx_user_symbols(i+3,1));conj(rx_user_symbols(i+3,2));
                rx_user_symbols(i,3);rx_user_symbols(i,4);conj(rx_user_symbols(i+1,3));conj(rx_user_symbols(i+1,4));
                rx_user_symbols(i+2,3);rx_user_symbols(i+2,4);conj(rx_user_symbols(i+3,3));conj(rx_user_symbols(i+3,4))];
            
            H_est_user_combined = [H_est_user(i,1,1),-H_est_user(i,1,3),0,0;H_est_user(i,2,1),-H_est_user(i,2,3),0,0;conj(H_est_user(i+1,1,3)),conj(H_est_user(i+1,1,1)),0,0;conj(H_est_user(i+1,2,3)),conj(H_est_user(i+1,2,1)),0,0;...  % generate equivalent channel matrix
                0,0,H_est_user(i+2,1,2),-H_est_user(i+2,1,4);0,0,H_est_user(i+2,2,2),-H_est_user(i+2,2,4);0,0,conj(H_est_user(i+3,1,4)),conj(H_est_user(i+3,1,2));0,0,conj(H_est_user(i+3,2,4)),conj(H_est_user(i+3,2,2));...
                H_est_user(i,3,1),-H_est_user(i,3,3),0,0;H_est_user(i,4,1),-H_est_user(i,4,3),0,0;conj(H_est_user(i+1,3,3)),conj(H_est_user(i+1,3,1)),0,0;conj(H_est_user(i+1,4,3)),conj(H_est_user(i+1,4,1)),0,0;...  % generate equivalent channel matrix
                0,0,H_est_user(i+2,3,2),-H_est_user(i+2,3,4);0,0,H_est_user(i+2,4,2),-H_est_user(i+2,4,4);0,0,conj(H_est_user(i+3,3,4)),conj(H_est_user(i+3,3,2));0,0,conj(H_est_user(i+3,4,4)),conj(H_est_user(i+3,4,2))];
            switch receiver
                case 'SSD'
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2(1:2) = pinv(H_est_user_combined([1:4,9:12],1:2))*rx_user_symbols_combined([1:4,9:12])*sqrt(2);    % calculate received symbols
                    rx_user_symbols_combined_2(3:4) = pinv(H_est_user_combined([5:8,13:16],3:4))*rx_user_symbols_combined([5:8,13:16])*sqrt(2);
%                     rx_user_symbols_combined_2(5:6) = pinv(H_est_user_combined(9:12,1:2))*rx_user_symbols_combined(9:12)*sqrt(2);    % calculate received symbols
%                     rx_user_symbols_combined_2(7:8) = pinv(H_est_user_combined(13:16,3:4))*rx_user_symbols_combined(13:16)*sqrt(2);
                    
                    rx_layer_x(:,(i+3)/4)=rx_user_symbols_combined_2;
                    [Q,R]= qr(H_est_user_combined([1:4,9:12],1:2));
                    LLR_SD(1:sum(M(1:2)),(i+3)/4) = LTE_softsphere(rx_layer_x(1:2,(i+3)/4),rx_user_symbols_combined([1:4,9:12]).'*sqrt(2),Q,R,symbol_alphabet(1:2,:),bittable(1:sum(M(1:2)),:),nLayers/2,M(1:2));
                    
                    [Q,R]= qr(H_est_user_combined([5:8,13:16],3:4));
                    LLR_SD(sum(M(1:2))+1:sum(M),(i+3)/4) = LTE_softsphere(rx_layer_x(3:4,(i+3)/4),rx_user_symbols_combined([5:8,13:16]).'*sqrt(2),Q,R,symbol_alphabet(3:4,:),bittable(1:sum(M(1:2)),:),nLayers/2,M(3:4));
                    
                case 'ZF'
                    inv_temp = pinv(H_est_user_combined);
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+3)/4)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+3)/4)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x(3,(i+3)/4)=rx_user_symbols_combined_2(3,:);
                    rx_layer_x(4,(i+3)/4)=rx_user_symbols_combined_2(4,:);
                    %%%%%%%%%%%%%%%%%%%
                    
                    rx_layer_x_equalized((i+3)/4,1) = rx_user_symbols_combined_2(1,:);
                    rx_layer_x_equalized((i+3)/4,2) = rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+3)/4,3) = rx_user_symbols_combined_2(3,:);
                    rx_layer_x_equalized((i+3)/4,4) = rx_user_symbols_combined_2(4,:)';
                    
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1);...
                        noise_enhancement_tmp(3)*ones(M(3),1);noise_enhancement_tmp(4)*ones(M(4),1)];
                    LLR_SD(:,(i+3)/4) = LTE_demapper(rx_layer_x(:,(i+3)/4),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                case 'MMSE'
                    temp = H_est_user_combined'*H_est_user_combined;
                    inv_temp = (temp+sigma_n2*eye(size(temp)))^-1*H_est_user_combined';
                    % MMSE Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+3)/4)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+3)/4)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x(3,(i+3)/4)=rx_user_symbols_combined_2(3,:);
                    rx_layer_x(4,(i+3)/4)=rx_user_symbols_combined_2(4,:);
                    %%%%%%%%%%%%%%%%%%%
                    
                    rx_layer_x_equalized((i+3)/4,1) = rx_user_symbols_combined_2(1,:);
                    rx_layer_x_equalized((i+3)/4,2) = rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+3)/4,3) = rx_user_symbols_combined_2(3,:);
                    rx_layer_x_equalized((i+3)/4,4) = rx_user_symbols_combined_2(4,:)';
                    
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1);...
                        noise_enhancement_tmp(3)*ones(M(3),1);noise_enhancement_tmp(4)*ones(M(4),1)];
                    LLR_SD(:,(i+3)/4) = LTE_demapper(rx_layer_x(:,(i+3)/4),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
            end
        end
    
    
    elseif (UE_nRX == 2)                % 4x2 case
        for i = 1:4:length(rx_user_symbols)-3
            rx_user_symbols_combined = [rx_user_symbols(i,1);rx_user_symbols(i,2);conj(rx_user_symbols(i+1,1));conj(rx_user_symbols(i+1,2));
                rx_user_symbols(i+2,1);rx_user_symbols(i+2,2);conj(rx_user_symbols(i+3,1));conj(rx_user_symbols(i+3,2))];
            H_est_user_combined = [H_est_user(i,1,1),-H_est_user(i,1,3),0,0;H_est_user(i,2,1),-H_est_user(i,2,3),0,0;conj(H_est_user(i+1,1,3)),conj(H_est_user(i+1,1,1)),0,0;conj(H_est_user(i+1,2,3)),conj(H_est_user(i+1,2,1)),0,0;...  % generate equivalent channel matrix
                0,0,H_est_user(i+2,1,2),-H_est_user(i+2,1,4);0,0,H_est_user(i+2,2,2),-H_est_user(i+2,2,4);0,0,conj(H_est_user(i+3,1,4)),conj(H_est_user(i+3,1,2));0,0,conj(H_est_user(i+3,2,4)),conj(H_est_user(i+3,2,2))];

            switch receiver
                case 'SSD'
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2(1:2) = pinv(H_est_user_combined(1:4,1:2))*rx_user_symbols_combined(1:4)*sqrt(2);    % calculate received symbols
                    rx_user_symbols_combined_2(3:4) = pinv(H_est_user_combined(5:8,3:4))*rx_user_symbols_combined(5:8)*sqrt(2);
                    
                    rx_layer_x(:,(i+3)/4)=rx_user_symbols_combined_2;
                    
                    [Q,R]= qr(H_est_user_combined(1:4,1:2));
                    LLR_SD(1:sum(M(1:2)),(i+3)/4) = LTE_softsphere(rx_layer_x(1:2,(i+3)/4),rx_user_symbols_combined(1:4).'*sqrt(2),Q,R,symbol_alphabet(1:2,:),bittable(1:sum(M(1:2)),:),nLayers/2,M(1:2));
                    
                    [Q,R]= qr(H_est_user_combined(5:8,3:4));
                    LLR_SD(sum(M(1:2))+1:sum(M),(i+3)/4) = LTE_softsphere(rx_layer_x(3:4,(i+3)/4),rx_user_symbols_combined(5:8).'*sqrt(2),Q,R,symbol_alphabet(3:4,:),bittable(1:sum(M(1:2)),:),nLayers/2,M(3:4));
                    
                case 'ZF'
                    inv_temp = pinv(H_est_user_combined);
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+3)/4)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+3)/4)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x(3,(i+3)/4)=rx_user_symbols_combined_2(3,:);
                    rx_layer_x(4,(i+3)/4)=rx_user_symbols_combined_2(4,:);
                    %%%%%%%%%%%%%%%%%%%
                    
                    rx_layer_x_equalized((i+3)/4,1) = rx_user_symbols_combined_2(1,:);
                    rx_layer_x_equalized((i+3)/4,2) = rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+3)/4,3) = rx_user_symbols_combined_2(3,:);
                    rx_layer_x_equalized((i+3)/4,4) = rx_user_symbols_combined_2(4,:)';
                    
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1);...
                        noise_enhancement_tmp(3)*ones(M(3),1);noise_enhancement_tmp(4)*ones(M(4),1)];
                    LLR_SD(:,(i+3)/4) = LTE_demapper(rx_layer_x(:,(i+3)/4),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                case 'MMSE'
                    temp = H_est_user_combined'*H_est_user_combined;
                    inv_temp = (temp+sigma_n2*eye(size(temp)))^-1*H_est_user_combined';
                    % MMSE Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+3)/4)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+3)/4)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x(3,(i+3)/4)=rx_user_symbols_combined_2(3,:);
                    rx_layer_x(4,(i+3)/4)=rx_user_symbols_combined_2(4,:);
                    %%%%%%%%%%%%%%%%%%%
                    
                    rx_layer_x_equalized((i+3)/4,1) = rx_user_symbols_combined_2(1,:);
                    rx_layer_x_equalized((i+3)/4,2) = rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+3)/4,3) = rx_user_symbols_combined_2(3,:);
                    rx_layer_x_equalized((i+3)/4,4) = rx_user_symbols_combined_2(4,:)';
                    
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1);...
                        noise_enhancement_tmp(3)*ones(M(3),1);noise_enhancement_tmp(4)*ones(M(4),1)];
                    LLR_SD(:,(i+3)/4) = LTE_demapper(rx_layer_x(:,(i+3)/4),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
            end
        end
        
    elseif (UE_nRX == 1)            % 4x1 case
        for i = 1:4:length(rx_user_symbols)-3
            rx_user_symbols_combined = [rx_user_symbols(i,1);conj(rx_user_symbols(i+1,1));rx_user_symbols(i+2,1);conj(rx_user_symbols(i+3,1))];  % generate equivalent channel matrix
           H_est_user_combined =    [H_est_user(i,1,1),-H_est_user(i,1,3),0,0;conj(H_est_user(i+1,1,3)),conj(H_est_user(i+1,1,1)),0,0;...
                0,0,H_est_user(i+2,1,2),-H_est_user(i+2,1,4);0,0,conj(H_est_user(i+3,1,4)),conj(H_est_user(i+3,1,2))];
            
            switch receiver
                case 'SSD'
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2(1:2) = pinv(H_est_user_combined(1:2,1:2))*rx_user_symbols_combined(1:2)*sqrt(2);    % calculate received symbols
                    rx_user_symbols_combined_2(3:4) = pinv(H_est_user_combined(3:4,3:4))*rx_user_symbols_combined(3:4)*sqrt(2);
                    
                    rx_layer_x(:,(i+3)/4)=rx_user_symbols_combined_2;
                    
                    %%%%%%%%%%%%%%%%%%%
                    [Q,R]= qr(H_est_user_combined(1:2,1:2));
                    LLR_SD(1:sum(M(1:2)),(i+3)/4) = LTE_softsphere(rx_layer_x(1:2,(i+3)/4),rx_user_symbols_combined(1:2).'*sqrt(2),Q,R,symbol_alphabet(1:2,:),bittable(1:sum(M(1:2)),:),nLayers/2,M(1:2));
                    
                    [Q,R]= qr(H_est_user_combined(3:4,3:4));
                    LLR_SD(sum(M(1:2))+1:sum(M),(i+3)/4) = LTE_softsphere(rx_layer_x(3:4,(i+3)/4),rx_user_symbols_combined(3:4).'*sqrt(2),Q,R,symbol_alphabet(3:4,:),bittable(1:sum(M(1:2)),:),nLayers/2,M(3:4));
                    
                case 'ZF'
                    inv_temp = pinv(H_est_user_combined);
                    % ZF Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+3)/4)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+3)/4)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x(3,(i+3)/4)=rx_user_symbols_combined_2(3,:);
                    rx_layer_x(4,(i+3)/4)=rx_user_symbols_combined_2(4,:);
                    %%%%%%%%%%%%%%%%%%%
                    
                    rx_layer_x_equalized((i+3)/4,1) = rx_user_symbols_combined_2(1,:);
                    rx_layer_x_equalized((i+3)/4,2) = rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+3)/4,3) = rx_user_symbols_combined_2(3,:);
                    rx_layer_x_equalized((i+3)/4,4) = rx_user_symbols_combined_2(4,:)';
                    
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1);...
                        noise_enhancement_tmp(3)*ones(M(3),1);noise_enhancement_tmp(4)*ones(M(4),1)];
                    LLR_SD(:,(i+3)/4) = LTE_demapper(rx_layer_x(:,(i+3)/4),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
                case 'MMSE'
                    temp = H_est_user_combined'*H_est_user_combined;
                    inv_temp = (temp+sigma_n2*eye(size(temp)))^-1*H_est_user_combined';
                    % MMSE Detector
                    %%%%%%%%%%%%%%%%%%%
                    rx_user_symbols_combined_2 = inv_temp*rx_user_symbols_combined*sqrt(2);    % calculate received symbols
                    rx_layer_x(1,(i+3)/4)=rx_user_symbols_combined_2(1,:);                       % map them to the corresponding layer
                    rx_layer_x(2,(i+3)/4)=rx_user_symbols_combined_2(2,:);
                    rx_layer_x(3,(i+3)/4)=rx_user_symbols_combined_2(3,:);
                    rx_layer_x(4,(i+3)/4)=rx_user_symbols_combined_2(4,:);
                    %%%%%%%%%%%%%%%%%%%
                    
                    rx_layer_x_equalized((i+3)/4,1) = rx_user_symbols_combined_2(1,:);
                    rx_layer_x_equalized((i+3)/4,2) = rx_user_symbols_combined_2(2,:)';
                    rx_layer_x_equalized((i+3)/4,3) = rx_user_symbols_combined_2(3,:);
                    rx_layer_x_equalized((i+3)/4,4) = rx_user_symbols_combined_2(4,:)';
                    
                    Hg = inv_temp*H_est_user_combined;
                    noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                    noise_enhancement = [noise_enhancement_tmp(1)*ones(M(1),1);noise_enhancement_tmp(2)*ones(M(2),1);...
                        noise_enhancement_tmp(3)*ones(M(3),1);noise_enhancement_tmp(4)*ones(M(4),1)];
                    LLR_SD(:,(i+3)/4) = LTE_demapper(rx_layer_x(:,(i+3)/4),symbol_alphabet,bittable,nLayers,M,Hg,noise_enhancement);
            end
        end
    else error('number of RX antennas not supported');
    end
end