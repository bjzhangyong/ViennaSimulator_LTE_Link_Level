function [channel_estimate] = LTE_channel_estimator(LTE_params,ChanMod, UE, rx_ref_symbols, RefSym, RefMapping, perfect_channel, sigma_n2, y_tx_assembled_genie, y_rx_assembled)
% LTE channel estimator - to filter the output of the transmitter.
% [chan_output] = LTE_channel_model(BS_output, SNR)
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) by INTHFT
% www.nt.tuwien.ac.at
%
% input :   ChanMod                     ... struct -> include info about filtering type, channel autocorrelation matrix
%           channel_estimation_method   ... string -> type of channel
%           estimation
%           channel_interpolation_method... string -> type of channel interpolation
%           rx_ref_symbols              ... [ x ] matrix of recieved reference symbols
%           RefSym                      ... [ x ] matrix of transmitted reference symbols
%           RefMapping                  ... [ x ] logical matrix of positions of reference symbols
%           perfect_channel             ... [ x ] coefficients of perfect channel
%           sigma_n2                    ... [1x1] double noise estimate
% output:   channel_estimate            ... [ x ] matrix of estimate of channel coefficients
%
% date of creation: 2009/04/06
% last changes: 2008/12/12  Simko

% global LTE_params;

%% Channel estimator
channel_estimate = nan(LTE_params.Ntot,LTE_params.Nsub,ChanMod.nRX,ChanMod.nTX);
interpolation = false;
channel_estimation_method = UE.channel_estimation_method;
channel_interpolation_method = UE.channel_interpolation_method;

switch ChanMod.filtering
    case 'BlockFading'
        switch channel_estimation_method
            case 'corrMMSE'
                [position_freq,position_time] = find(RefMapping(:,:,1));
                pilots_num = length(position_freq);
                H_ls_matrix = nan(pilots_num/2,ChanMod.nRX,ChanMod.nTX);
                for tt = 1:ChanMod.nTX
                    for rr = 1:ChanMod.nRX
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq),2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq)/2,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_matrix(:,rr,tt) = H_ls_one_channel(reference_mapping);
                    end
                end
                H_ls_vector = H_ls_matrix(:);
                channel_autocorrelation = UE.channel_autocorrelation_matrix;
                antenna_correlation = kron(squeeze(ChanMod.corrRX(1,:,:)) , squeeze(ChanMod.corrTX(1,:,:)));
                channel_autocorrelation_full = kron(antenna_correlation,channel_autocorrelation);
                permutation_matrix = zeros(LTE_params.Ntot*ChanMod.nTX*ChanMod.nRX);
                permutation_vector = 1:LTE_params.Ntot*ChanMod.nTX*ChanMod.nRX;
                pilot_pos = repmat(notused,ChanMod.nRX*ChanMod.nTX,1) + kron((0:(ChanMod.nRX*ChanMod.nTX-1))',ones(LTE_params.Ntot/3,1)*LTE_params.Ntot);
                pilots_num = length(pilot_pos);
                permutation_vector(pilot_pos) = [];
                permutation_vector = [pilot_pos;permutation_vector'];
                for pos_i = 1:LTE_params.Ntot*ChanMod.nRX*ChanMod.nTX
                    permutation_matrix(pos_i,permutation_vector(pos_i)) = 1;
                end
                channel_autocorrelation_full_permutated = permutation_matrix * channel_autocorrelation_full * permutation_matrix';
                channel_estimate_permutated = channel_autocorrelation_full_permutated(:,1:pilots_num)*inv( channel_autocorrelation_full_permutated(1:pilots_num,1:pilots_num) + sigma_n2*eye(pilots_num)) * H_ls_vector;
                channel_estimate_help = permutation_matrix'*channel_estimate_permutated;
                channel_estimate_help2 = reshape(channel_estimate_help,LTE_params.Ntot,ChanMod.nRX,ChanMod.nTX);
                channel_estimate(:,1,:,:) = channel_estimate_help2;
                channel_estimate = repmat(channel_estimate(:,1,:,:),[1,LTE_params.Nsub,1,1]);
        end
end

for tt = 1:ChanMod.nTX
    for rr = 1:ChanMod.nRX
        switch ChanMod.filtering
            case 'BlockFading'
                switch channel_estimation_method
                    case 'PERFECT'
                        channel_estimate = perfect_channel; % perfect channel knowledge
                    case 'LS' %LS for block fading
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq),2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq)/2,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_one_channel = H_ls_one_channel(reference_mapping);
                        interpolation = true;
                    case 'MMSE'
                        %first LS estimate on teh pilots positions
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            pilots_num = length(position_freq);
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,pilots_num,2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            pilots_num = length(position_freq)/2;
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,pilots_num,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_one_channel = H_ls_one_channel(reference_mapping);
                        
                        %MMSE estimate in case of ideal autocorrelation matrix or
                        %the estimation of the autocorrelation matrix is ready
                        if(strcmp(UE.autocorrelation_matrix_type,'ideal') || UE.realization_num >= UE.realization_num_total)
                            channel_autocorrelation = UE.channel_autocorrelation_matrix;
                            permutation_matrix = zeros(LTE_params.Ntot);
                            permutation_vector = 1:LTE_params.Ntot;
                            permutation_vector(position(:,1)) = [];
                            permutation_vector = [sort(position(:,1));permutation_vector'];
                            [~,permutation_vector_back] = sort(permutation_vector);
                            channel_autocorrelation_permutated = channel_autocorrelation(permutation_vector,permutation_vector);
                            channel_estimate_permutated = channel_autocorrelation_permutated(1:LTE_params.Ntot,1:pilots_num)*inv( channel_autocorrelation_permutated(1:pilots_num,1:pilots_num) + sigma_n2*eye(pilots_num)) * H_ls_one_channel;
                            channel_estimate_help = channel_estimate_permutated(permutation_vector_back);
                        end
                        %estimation of autocorrelation matrix
                        if( strcmp(UE.autocorrelation_matrix_type,'estimated') && UE.realization_num < UE.realization_num_total )
                            %linear interpolation
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            channel_estimate_help = interp1(first_ref_sym:3:LTE_params.Ntot,H_ls_one_channel,1:LTE_params.Ntot,'linear','extrap');
                            channel_estimate_help = channel_estimate_help.';
                            %update of autocorrelation amtrix
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix * UE.realization_num + channel_estimate_help*channel_estimate_help';
                            UE.realization_num = UE.realization_num + 1;
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix/UE.realization_num;
                        end
                        channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help,1 , LTE_params.Nsub);
                    case 'ALMMSE'
                        %Mehlfuhrer, C.; Caban, S.; Rupp, M., "An accurate and low complex channel estimator for OFDM WiMAX," Communications, Control and Signal Processing, 2008. ISCCSP 2008. 3rd International Symposium on , vol., no., pp.922-926, 12-14 March 2008
                        %URL:http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=45
                        %37355&isnumber=4537177
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq),2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq)/2,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_one_channel = H_ls_one_channel(reference_mapping);
                        channel_estimate_help = nan(LTE_params.Ntot,1);
                        
                        length_2 = 12;   %number of adjacent subcarriers, which are used for channel estimation
                        channel_autocorrelation = UE.channel_autocorrelation_matrix;
                        channel_autocorrelation_small = zeros(length_2);
                        for step_i=1:LTE_params.Ntot/length_2
                            channel_autocorrelation_small = channel_autocorrelation_small + channel_autocorrelation((step_i-1)*length_2 + 1:(step_i)*length_2,(step_i-1)*length_2 + 1:(step_i)*length_2);
                        end
                        channel_autocorrelation_small = channel_autocorrelation_small/(LTE_params.Ntot/length_2);
                        
                        %estimte for the first position
                        permutation_matrix = zeros(length_2);
                        permutation_vector = 1:length_2;
                        pilots_small = notused <= length_2;
                        permutation_vector(notused(pilots_small)) = [];
                        permutation_vector = [notused(pilots_small);permutation_vector'];
                        for pos_i = 1:length_2
                            permutation_matrix(pos_i,permutation_vector(pos_i)) = 1;
                        end
                        pilots_small_num = length(notused(pilots_small));
                        channel_autocorrelation_small_permutated = permutation_matrix * channel_autocorrelation_small * permutation_matrix';
                        estimate_small_permuted = channel_autocorrelation_small_permutated(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num)) * H_ls_one_channel(pilots_small);
                        estimate_small = permutation_matrix'*estimate_small_permuted;
                        channel_estimate_help(1:floor((length_2+1)/2)) = estimate_small(1:floor((length_2+1)/2));
                        %estimate for the middle position
                        for channelpos_i=floor((length_2+1)/2) + 1:LTE_params.Ntot-floor((length_2 -1)/2)-1
                            %for channelpos_i=floor((length_2+1)/2) + 1:LTE_params.Ntot-floor((length_2 -1)/2)
                            permutation_matrix = zeros(length_2);
                            permutation_vector = 1:length_2;
                            start = channelpos_i - fix((length_2-1)/2);
                            stop = channelpos_i + ceil((length_2-1)/2);
                            pilots_small = start <= notused & notused <= stop;
                            pilots_pos_small = notused(pilots_small) - start + 1;
                            %pilots_pos_small = mod(notused(pilots_small),length_2);
                            %pilots_pos_small(find(pilots_pos_small == 0)) = pilots_pos_small(find(pilots_pos_small == 0)) + length_2;
                            permutation_vector(pilots_pos_small) = [];
                            permutation_vector = [pilots_pos_small;permutation_vector'];
                            for pos_i = 1:length_2
                                permutation_matrix(pos_i,permutation_vector(pos_i)) = 1;
                            end
                            pilots_small_num = length(pilots_pos_small);
                            channel_autocorrelation_small_permutated = permutation_matrix * channel_autocorrelation_small * permutation_matrix';
                            F = channel_autocorrelation_small_permutated(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num));
                            estimate_small_permuted = F * H_ls_one_channel(pilots_small);
                            estimate_small = permutation_matrix'*estimate_small_permuted;
                            %channel_estimate(channelpos_i) = estimate_small(ceil((length_2+1)/2));
                            channel_estimate_help(channelpos_i) = estimate_small(channelpos_i - start +1 );
                        end
                        %estimate for last position
                        permutation_matrix = zeros(length_2);
                        permutation_vector = 1:length_2;
                        pilots_small = notused >= LTE_params.Ntot-length_2+1;
                        
                        pilots_pos_small = notused(pilots_small) - (LTE_params.Ntot-length_2);
                        if(max(pilots_pos_small)>length_2)
                            pilots_pos_small = mod(notused(pilots_small),length_2);
                            pilots_pos_small(find(pilots_pos_small == 0)) = pilots_pos_small(find(pilots_pos_small == 0)) + length_2;
                        end
                        permutation_vector(pilots_pos_small) = [];
                        permutation_vector = [pilots_pos_small;permutation_vector'];
                        for pos_i = 1:length_2
                            permutation_matrix(pos_i,permutation_vector(pos_i)) = 1;
                        end
                        pilots_small_num = length(notused(pilots_small));
                        channel_autocorrelation_small_permutated = permutation_matrix * channel_autocorrelation_small * permutation_matrix';
                        estimate_small_permuted = channel_autocorrelation_small_permutated(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num)) * H_ls_one_channel(pilots_small);
                        estimate_small = permutation_matrix'*estimate_small_permuted;
                        channel_estimate_help(end-fix((length_2-1)/2):end) = estimate_small(fix((length_2+1)/2)+1:end);
                        
                        channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help,1 , LTE_params.Nsub);
                        interpolation = false;
                    case 'corrMMSE'
                        
                    case 'ALMMSE2'
                        %Mehlfuhrer, C.; Caban, S.; Rupp, M., "An accurate and low complex channel estimator for OFDM WiMAX," Communications, Control and Signal Processing, 2008. ISCCSP 2008. 3rd International Symposium on , vol., no., pp.922-926, 12-14 March 2008
                        %URL:http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=45
                        %37355&isnumber=4537177
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq),2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq)/2,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_one_channel = H_ls_one_channel(reference_mapping);
                        channel_estimate_help = nan(LTE_params.Ntot,1);
                        
                        if( strcmp(UE.autocorrelation_matrix_type,'estimated') && UE.realization_num < UE.realization_num_total )
                            %linear interpolation
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            channel_estimate_help = interp1(first_ref_sym:3:LTE_params.Ntot,H_ls_one_channel,1:LTE_params.Ntot,'linear','extrap');
                            channel_estimate_help = channel_estimate_help.';
                            %update of autocorrelation amtrix
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix * UE.realization_num + channel_estimate_help*channel_estimate_help';
                            UE.realization_num = UE.realization_num + 1;
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix/UE.realization_num;
                        end
                        
                        if(strcmp(UE.autocorrelation_matrix_type,'ideal') || UE.realization_num >= UE.realization_num_total)
                            length_2 =12;   %number of adjacent subcarriers, which are used for channel estimation
                            channel_autocorrelation = UE.channel_autocorrelation_matrix;
                            channel_autocorrelation_small = zeros(length_2);
                            length_shift = floor(length_2/2);
                            for step_i=1:LTE_params.Ntot/length_shift-1
                                channel_autocorrelation_small = channel_autocorrelation_small + channel_autocorrelation((step_i-1)*length_shift + 1:(step_i-1)*length_shift + length_2,(step_i-1)*length_shift + 1:(step_i-1)*length_shift + length_2);
                            end
                            channel_autocorrelation_small = channel_autocorrelation_small/(LTE_params.Ntot/length_shift-1);
                            
                            %calculation of three different filtering matrices
                            permutation_vector1 = 1:length_2;
                            permutation_matrix1 = zeros(length_2);
                            permutation_vector1(1:3:length_2) = [];
                            permutation_vector1 = [1:3:length_2,permutation_vector1];
                            for pos_i = 1:length_2
                                permutation_matrix1(pos_i,permutation_vector1(pos_i)) = 1;
                            end
                            
                            permutation_vector2 = 1:length_2;
                            permutation_matrix2 = zeros(length_2);
                            permutation_vector2(2:3:length_2) = [];
                            permutation_vector2 = [2:3:length_2,permutation_vector2];
                            for pos_i = 1:length_2
                                permutation_matrix2(pos_i,permutation_vector2(pos_i)) = 1;
                            end
                            
                            permutation_vector3 = 1:length_2;
                            permutation_matrix3 = zeros(length_2);
                            permutation_vector3(3:3:length_2) = [];
                            permutation_vector3 = [3:3:length_2,permutation_vector3];
                            for pos_i = 1:length_2
                                permutation_matrix3(pos_i,permutation_vector3(pos_i)) = 1;
                            end
                            
                            pilots_small_num = floor(length_2/3);
                            
                            channel_autocorrelation_small_permutated1 = permutation_matrix1 * channel_autocorrelation_small * permutation_matrix1';
                            F1 = channel_autocorrelation_small_permutated1(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated1(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num));
                            
                            channel_autocorrelation_small_permutated2 = permutation_matrix2 * channel_autocorrelation_small * permutation_matrix2';
                            F2 = channel_autocorrelation_small_permutated2(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated2(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num));
                            
                            channel_autocorrelation_small_permutated3 = permutation_matrix3 * channel_autocorrelation_small * permutation_matrix3';
                            F3 = channel_autocorrelation_small_permutated3(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated3(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num));
                            
                            %estimte for the first position
                            
                            pilots_small = notused <= length_2;
                            switch notused(1)
                                case 1
                                    F = F1;
                                    permutation_matrix = permutation_matrix1;
                                case 2
                                    F = F2;
                                    permutation_matrix = permutation_matrix2;
                                case 3
                                    F = F3;
                                    permutation_matrix = permutation_matrix3;
                            end
                            estimate_small_permuted = F * H_ls_one_channel(pilots_small);
                            estimate_small = permutation_matrix'*estimate_small_permuted;
                            channel_estimate_help(1:floor((length_2+1)/2)) = estimate_small(1:floor((length_2+1)/2));
                            
                            %estimate for the middle position
                            
                            for channelpos_i=floor((length_2+1)/2) + 1:LTE_params.Ntot-floor((length_2 -1)/2)-1
                                start = channelpos_i - fix((length_2-1)/2);
                                stop = channelpos_i + ceil((length_2-1)/2);
                                pilots_small = start <= notused & notused <= stop;
                                pilots_pos_small = notused(pilots_small) - start + 1;
                                
                                switch pilots_pos_small(1)
                                    case 1
                                        F = F1;
                                        permutation_matrix = permutation_matrix1;
                                    case 2
                                        F = F2;
                                        permutation_matrix = permutation_matrix2;
                                    case 3
                                        F = F3;
                                        permutation_matrix = permutation_matrix3;
                                end
                                
                                estimate_small_permuted = F * H_ls_one_channel(pilots_small);
                                estimate_small = permutation_matrix'*estimate_small_permuted;
                                %channel_estimate_help(channelpos_i) = estimate_small(channelpos_i - start +1 );
                                channel_estimate_help(channelpos_i) = estimate_small(floor((length_2+1)/2));
                            end
                            %estimate for last position
                            pilots_small = notused > LTE_params.Ntot - length_2;
                            pilots_pos_small = notused(pilots_small) - (LTE_params.Ntot - length_2);
                            switch pilots_pos_small(1)
                                case 1
                                    F = F1;
                                    permutation_matrix = permutation_matrix1;
                                case 2
                                    F = F2;
                                    permutation_matrix = permutation_matrix2;
                                case 3
                                    F = F3;
                                    permutation_matrix = permutation_matrix3;
                            end
                            estimate_small_permuted = F * H_ls_one_channel(pilots_small);
                            estimate_small = permutation_matrix'*estimate_small_permuted;
                            channel_estimate_help(end-fix((length_2-1)/2):end) = estimate_small(fix((length_2+1)/2)+1:end);
                        end
                        
                        channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help,1 , LTE_params.Nsub);
                        interpolation = false;
                    case 'ALMMSE3'
                        %Mehlfuhrer, C.; Caban, S.; Rupp, M., "An accurate and low complex channel estimator for OFDM WiMAX," Communications, Control and Signal Processing, 2008. ISCCSP 2008. 3rd International Symposium on , vol., no., pp.922-926, 12-14 March 2008
                        %URL:http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=45
                        %37355&isnumber=4537177
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq),2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq)/2,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_one_channel = H_ls_one_channel(reference_mapping);
                        channel_estimate_help = nan(LTE_params.Ntot,1);
                        
                        first_ref_sym = min(position_freq); %position of first refernce symbol
                        channel_estimate_help1 = interp1(first_ref_sym:3:LTE_params.Ntot,H_ls_one_channel,1:LTE_params.Ntot,channel_interpolation_method,'extrap');
                        channel_estimate_help1 = channel_estimate_help1.';
                        
                        length_2 =12;   %number of adjacent subcarriers, which are used for channel estimation
                        
                        %channel_autocorrelation = UE.channel_autocorrelation_matrix;
                        channel_autocorrelation = channel_estimate_help1*channel_estimate_help1';
                        channel_autocorrelation_small = zeros(length_2);
                        length_shift = floor(length_2/2);
                        for step_i=1:LTE_params.Ntot/length_shift-1
                            channel_autocorrelation_small = channel_autocorrelation_small + channel_autocorrelation((step_i-1)*length_shift + 1:(step_i-1)*length_shift + length_2,(step_i-1)*length_shift + 1:(step_i-1)*length_shift + length_2);
                        end
                        channel_autocorrelation_small = channel_autocorrelation_small/(LTE_params.Ntot/length_shift-1);
                        
                        %calculation of three different filtering matrices
                        permutation_vector1 = 1:length_2;
                        permutation_matrix1 = zeros(length_2);
                        permutation_vector1(1:3:length_2) = [];
                        permutation_vector1 = [1:3:length_2,permutation_vector1];
                        for pos_i = 1:length_2
                            permutation_matrix1(pos_i,permutation_vector1(pos_i)) = 1;
                        end
                        
                        permutation_vector2 = 1:length_2;
                        permutation_matrix2 = zeros(length_2);
                        permutation_vector2(2:3:length_2) = [];
                        permutation_vector2 = [2:3:length_2,permutation_vector2];
                        for pos_i = 1:length_2
                            permutation_matrix2(pos_i,permutation_vector2(pos_i)) = 1;
                        end
                        
                        permutation_vector3 = 1:length_2;
                        permutation_matrix3 = zeros(length_2);
                        permutation_vector3(3:3:length_2) = [];
                        permutation_vector3 = [3:3:length_2,permutation_vector3];
                        for pos_i = 1:length_2
                            permutation_matrix3(pos_i,permutation_vector3(pos_i)) = 1;
                        end
                        
                        pilots_small_num = floor(length_2/3);
                        
                        channel_autocorrelation_small_permutated1 = permutation_matrix1 * channel_autocorrelation_small * permutation_matrix1';
                        F1 = channel_autocorrelation_small_permutated1(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated1(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num));
                        
                        channel_autocorrelation_small_permutated2 = permutation_matrix2 * channel_autocorrelation_small * permutation_matrix2';
                        F2 = channel_autocorrelation_small_permutated2(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated2(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num));
                        
                        channel_autocorrelation_small_permutated3 = permutation_matrix3 * channel_autocorrelation_small * permutation_matrix3';
                        F3 = channel_autocorrelation_small_permutated3(:,1:pilots_small_num)*inv(  channel_autocorrelation_small_permutated3(1:pilots_small_num,1:pilots_small_num) + sigma_n2*eye(pilots_small_num));
                        
                        %estimte for the first position
                        
                        pilots_small = notused <= length_2;
                        switch notused(1)
                            case 1
                                F = F1;
                                permutation_matrix = permutation_matrix1;
                            case 2
                                F = F2;
                                permutation_matrix = permutation_matrix2;
                            case 3
                                F = F3;
                                permutation_matrix = permutation_matrix3;
                        end
                        estimate_small_permuted = F * H_ls_one_channel(pilots_small);
                        estimate_small = permutation_matrix'*estimate_small_permuted;
                        channel_estimate_help(1:floor((length_2+1)/2)) = estimate_small(1:floor((length_2+1)/2));
                        
                        %estimate for the middle position
                        
                        for channelpos_i=floor((length_2+1)/2) + 1:LTE_params.Ntot-floor((length_2 -1)/2)-1
                            start = channelpos_i - fix((length_2-1)/2);
                            stop = channelpos_i + ceil((length_2-1)/2);
                            pilots_small = start <= notused & notused <= stop;
                            pilots_pos_small = notused(pilots_small) - start + 1;
                            
                            switch pilots_pos_small(1)
                                case 1
                                    F = F1;
                                    permutation_matrix = permutation_matrix1;
                                case 2
                                    F = F2;
                                    permutation_matrix = permutation_matrix2;
                                case 3
                                    F = F3;
                                    permutation_matrix = permutation_matrix3;
                            end
                            
                            estimate_small_permuted = F * H_ls_one_channel(pilots_small);
                            estimate_small = permutation_matrix'*estimate_small_permuted;
                            %channel_estimate_help(channelpos_i) = estimate_small(channelpos_i - start +1 );
                            channel_estimate_help(channelpos_i) = estimate_small(floor((length_2+1)/2));
                        end
                        %estimate for last position
                        pilots_small = notused > LTE_params.Ntot - length_2;
                        pilots_pos_small = notused(pilots_small) - (LTE_params.Ntot - length_2);
                        switch pilots_pos_small(1)
                            case 1
                                F = F1;
                                permutation_matrix = permutation_matrix1;
                            case 2
                                F = F2;
                                permutation_matrix = permutation_matrix2;
                            case 3
                                F = F3;
                                permutation_matrix = permutation_matrix3;
                        end
                        estimate_small_permuted = F * H_ls_one_channel(pilots_small);
                        estimate_small = permutation_matrix'*estimate_small_permuted;
                        channel_estimate_help(end-fix((length_2-1)/2):end) = estimate_small(fix((length_2+1)/2)+1:end);
                        
                        channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help,1 , LTE_params.Nsub);
                        interpolation = false;
                    case  {'GENIE','GENIE_ALMMSE'}
                        %do nothing
                    otherwise
                        error('channel estimation not supported');
                end
                
                if(interpolation)
                    switch channel_interpolation_method
                        case {'linear', 'cubic', 'spline'}
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            channel_estimate_help = interp1(first_ref_sym:3:LTE_params.Ntot,H_ls_one_channel,1:LTE_params.Ntot,channel_interpolation_method,'extrap');
                            channel_estimate_help = channel_estimate_help.';
                        case 'sinc_freq'
                            pulse = zeros(LTE_params.Ntot,1);
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            pulse(first_ref_sym:3:end) = H_ls_one_channel;
                            channel_estimate_help = conv(pulse,sinc((-LTE_params.Ntot:LTE_params.Ntot)/3));
                            channel_estimate_help = channel_estimate_help((LTE_params.Ntot +1 ):(LTE_params.Ntot + LTE_params.Ntot)).';
                        case 'sinc_time'
                            pulse = zeros(LTE_params.Ntot,1);
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            pulse(first_ref_sym:3:end) = H_ls_one_channel;
                            PULSE = ifft(pulse,LTE_params.Nfft);
                            N = 3*[ones(floor(LTE_params.Nfft/(3*2)),1);zeros(LTE_params.Nfft - 2*floor(LTE_params.Nfft/(3*2)),1);ones(floor(LTE_params.Nfft/(3*2)),1)]; %rect
                            channel_time = PULSE.*N;
                            channel_estimate_help = fft(channel_time,LTE_params.Nfft);
                            channel_estimate_help = channel_estimate_help(1:LTE_params.Ntot);
                        case 'T-F'
                            first_ref_sym = min(position_freq); %position of first refernce symbol
                            channel_estimate_help = interpft(H_ls_one_channel,LTE_params.Ntot);
                            channel_estimate_help = circshift(channel_estimate_help,first_ref_sym-1);
                            
                        otherwise
                            error('interpolation method not supported')
                    end
                    channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help,1 , LTE_params.Nsub);
                end
            case 'FastFading'
                switch channel_estimation_method
                    case 'PERFECT'
                        channel_estimate = perfect_channel; % perfect channel knowledge
                    case 'LS'
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        interpolation = true;
                    case 'LS_block'
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        if(tt>2)
                            %position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position_freq); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq),2);
                            H_ls_one_channel(:,2) = [];
                        else
                            position = reshape(position_freq,length(position_freq)/2,2); %form the position of reference symbols to appropriate matrix
                            [notused reference_mapping] = sort(position(:,1)); %find the real order of reference symbols in frequency
                            H_ls_one_channel = reshape(H_ls_one_channel,length(position_freq)/2,2);
                        end
                        H_ls_one_channel = mean(H_ls_one_channel,2);
                        H_ls_one_channel = H_ls_one_channel(reference_mapping);
                        %interpolation
                        first_ref_sym = min(position_freq); %position of first refernce symbol
                        channel_estimate_help = interp1(first_ref_sym:3:LTE_params.Ntot,H_ls_one_channel,1:LTE_params.Ntot,channel_interpolation_method,'extrap');
                        channel_estimate_help = channel_estimate_help.';
                        channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help,1 , LTE_params.Nsub);
                        interpolation = false;
                    case 'MMSE_Rayleigh2'
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        if(tt>2)
                            H_ls_one_channel(:,3:4) = [];
                        end
                        [position_freq,position_time] = find(RefMapping(:,:,tt));  %positons of reference symbols
                        position_time_reshape = reshape(position_time,LTE_params.Ntot/6,[]);
                        position_time_small = position_time_reshape(1,:);
                        H_ls_mean_over_freq = mean(H_ls_one_channel,1);
                        
                        c = LTE_params.speed_of_light;
                        f = LTE_params.carrier_freq;  % Frequency at which our system operates
                        v = UE.user_speed;  %speed at which we move
                        w_d = 2*pi*v*f/c;   % maximum radian Doppler frequency
                        time_autocorrelation = toeplitz(besselj(0,w_d*(0:13)*137*LTE_params.SamplingTime));
                        
                        permutation_matrix = zeros(LTE_params.Nsub);
                        permutation_vector = 1:LTE_params.Nsub;
                        permutation_vector(position_time_small) = [];
                        permutation_vector = [position_time_small,permutation_vector];
                        for pos_i = 1:LTE_params.Nsub
                            permutation_matrix(pos_i,permutation_vector(pos_i)) = 1;
                        end
                        pilots_num = length(position_time_small);
                        channel_autocorrelation_permutated = permutation_matrix * time_autocorrelation * permutation_matrix';
                        channel_estimate_help_permutated = channel_autocorrelation_permutated(:,1:pilots_num)*inv( channel_autocorrelation_permutated(1:pilots_num,1:pilots_num) + sigma_n2*eye(pilots_num)) * H_ls_mean_over_freq.';
                        channel_estimate_help = permutation_matrix'*channel_estimate_help_permutated;
                        
                        channel_estimate(:,:,rr,tt) = repmat(channel_estimate_help.',[LTE_params.Ntot,1]);
                        
                        interpolation = false;
                    case 'MMSE'
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        H_ls_one_channel_vector = H_ls_one_channel(:);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));
                        if(tt>2)
                            pilots_num = length(position_freq);
                            position_freq_matrix = reshape(position_freq,pilots_num/2,2);
                            position_time_matrix = reshape(position_time,pilots_num/2,2);
                            H_ls_one_channel_vector(pilots_num + 1:end) = [];
                        else
                            pilots_num = length(position_freq);
                            position_freq_matrix = reshape(position_freq,pilots_num/4,4);
                            position_time_matrix = reshape(position_time,pilots_num/4,4);
                        end
                        
                        
                        c = LTE_params.speed_of_light;
                        f = LTE_params.carrier_freq;  % Frequency at which our system operates
                        v = UE.user_speed;  %speed at which we move
                        w_d = 2*pi*v*f/c;   % maximum radian Doppler frequency
                        time_autocorrelation = toeplitz(besselj(0,w_d*(0:13)*137*LTE_params.SamplingTime));
                        
                        if(strcmp(UE.autocorrelation_matrix_type,'ideal') || UE.realization_num >= UE.realization_num_total)
                            channel_autocorrelation = UE.channel_autocorrelation_matrix;
                            channel_autocorrelation_full = kron(time_autocorrelation,channel_autocorrelation);
                            permutation_matrix = zeros(LTE_params.Ntot*LTE_params.Nsub);
                            permutation_vector = 1:(LTE_params.Ntot*LTE_params.Nsub);
                            position = position_freq_matrix + LTE_params.Ntot*(position_time_matrix-1);
                            position = position(:);
                            permutation_vector(position) = [];
                            permutation_vector = [position;permutation_vector'];
                            [~,permutation_vector_back] = sort(permutation_vector);
                            channel_autocorrelation_permutated = channel_autocorrelation_full(permutation_vector,permutation_vector);
                            channel_estimate_help_permutated = channel_autocorrelation_permutated(:,1:pilots_num)*inv( channel_autocorrelation_permutated(1:pilots_num,1:pilots_num) + sigma_n2*eye(pilots_num)) * H_ls_one_channel_vector;
                            channel_estimate_help = channel_estimate_help_permutated(permutation_vector_back);
                            channel_estimate(:,:,rr,tt) = reshape(channel_estimate_help,LTE_params.Ntot,LTE_params.Nsub);
                        end
                        
                        if( strcmp(UE.autocorrelation_matrix_type,'estimated') && UE.realization_num < UE.realization_num_total )
                            %cubic interpolation
                            %add "virtual" pilots position, which we will extrapolate
                            if(tt>2)
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1), position_freq_matrix_ext(:,end), position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3, position_time_matrix_ext(:,end) + 6, position_time_matrix_ext(:,end) + 9];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            else
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            end
                            
                            %remove "virtual pilots" which are not nearest to real symbols
                            position_time_ext(position_not_needed) = [];
                            position_freq_ext(position_not_needed) = [];
                            position_to_extrapolate = zeros(length(position_time_ext)-length(position_time),2);
                            extrapolated_values = zeros(length(position_time_ext)-length(position_time),1);
                            three_nearest_points = zeros(length(position_time_ext)-length(position_time),3);
                            i = 1;
                            for pos_i = 1:length(position_time_ext)
                                if(sum(position_freq == position_freq_ext(pos_i) & position_time == position_time_ext(pos_i)))
                                else
                                    % here we find for every "virtual pilot" three
                                    % nearest pilots values and set the value of
                                    % "virtual pilot" as if it would be in plane,
                                    % which is spanned by those nearest points
                                    position_to_extrapolate(i,:) = [position_freq_ext(pos_i),position_time_ext(pos_i)];
                                    [distance , nearest_index] = sort((position_freq - position_freq_ext(pos_i)).^2 + (position_time - position_time_ext(pos_i)).^2);
                                    three_nearest_points(i,:) = nearest_index(1:3);
                                    three_nearest_points_help = sort(three_nearest_points(i,:));
                                    if(three_nearest_points_help == [three_nearest_points_help(2)-1,three_nearest_points_help(2),three_nearest_points_help(2)+1])
                                        positions_help = find(nearest_index < three_nearest_points_help(2) + 6 & nearest_index > three_nearest_points_help(2) - 6);
                                        nearest_index(positions_help) = [];
                                        three_nearest_points(i,3) = nearest_index(1);
                                    end
                                    grad = linsolve([position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,2)) , position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,3)) ; position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,2)) , position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,3))] ,[position_to_extrapolate(i,:)].' - [position_freq(three_nearest_points(i,1)); position_time(three_nearest_points(i,1))]);
                                    extrapolated_values(i) = H_ls_one_channel_vector(three_nearest_points(i,1)) + [H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,2)), H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,3))] * grad;
                                    i = i + 1;
                                end
                            end
                            [XI, YI] = meshgrid(1:LTE_params.Ntot,1:LTE_params.Nsub);
                            channel_estimate_help = griddata([position_freq; position_to_extrapolate(:,1)],[position_time; position_to_extrapolate(:,2)],[H_ls_one_channel_vector; extrapolated_values], XI, YI,'cubic');
                            channel_estimate(:,:,rr,tt) = channel_estimate_help.';
                            %update of autocorrelation amtrix
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix * UE.realization_num + (channel_estimate(:,:,rr,tt)*channel_estimate(:,:,rr,tt)')/LTE_params.Nsub;
                            UE.realization_num = UE.realization_num + 1;
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix/UE.realization_num;
                            %                             if UE.realization_num == UE.realization_num_total
                            %                                 correlation_ceof_vector_help = nan(LTE_params.Ntot,1);
                            %                                 for corr_i = 1:LTE_params.Ntot
                            %                                     matrix_help = UE.channel_autocorrelation_matrix(corr_i:LTE_params.Ntot,1:LTE_params.Ntot-corr_i+1);
                            %                                     correlation_ceof_vector_help(corr_i) = mean(diag(matrix_help));
                            %                                 end
                            %                                 UE.channel_autocorrelation_matrix = toeplitz(correlation_ceof_vector_help);
                            %                             end
                        end
                        interpolation = false;
                    case 'ALMMSE'
                        %M. Simko, C. Mehlführer, M. Wrulich, M. Rupp:"Doubly Dispersive Channel Estimation with Scalable Complexity";
                        %in: "Proceedings 2010 International ITG Workshop on Smart Antennas (WSA 2010)", (2010).
                        %URL:http://publik.tuwien.ac.at/files/PubDat_181229.pdf
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        H_ls_one_channel_vector = H_ls_one_channel(:);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));
                        if(tt>2)
                            pilots_num = length(position_freq);
                            position_freq_matrix = reshape(position_freq,pilots_num/2,2);
                            position_time_matrix = reshape(position_time,pilots_num/2,2);
                            H_ls_one_channel_vector(pilots_num + 1:end) = [];
                        else
                            pilots_num = length(position_freq);
                            position_freq_matrix = reshape(position_freq,pilots_num/4,4);
                            position_time_matrix = reshape(position_time,pilots_num/4,4);
                        end
                        
                        
                        c = LTE_params.speed_of_light;
                        f = LTE_params.carrier_freq;  % Frequency at which our system operates
                        v = UE.user_speed;  %speed at which we move
                        w_d = 2*pi*v*f/c;   % maximum radian Doppler frequency
                        time_autocorrelation = toeplitz(besselj(0,w_d*(0:13)*137*LTE_params.SamplingTime));
                        
                        if(strcmp(UE.autocorrelation_matrix_type,'ideal') || UE.realization_num >= UE.realization_num_total)
                            channel_autocorrelation = UE.channel_autocorrelation_matrix;
                            
                            
                            [eigenvector_t eigenvalue_t] = eig(time_autocorrelation);
                            eigenvalue_t = diag(eigenvalue_t);
                            [eigenvector_h eigenvalue_h] = eig(channel_autocorrelation);
                            eigenvalue_h = diag(eigenvalue_h);
                            
                            %removing small eigenvalues
                            treshold_time = LTE_params.eig_tresh_time;
                            treshold_freq = LTE_params.eig_tresh_freq;
                            
                            pos_t_small = find(abs(eigenvalue_t) < treshold_time);
                            pos_t_large = find(abs(eigenvalue_t) > treshold_time);
                            pos_h_small = find(abs(eigenvalue_h) < treshold_freq);
                            pos_h_large = find(abs(eigenvalue_h) > treshold_freq);
                            eigenvalue_t(pos_t_small) = [];
                            eigenvalue_h(pos_h_small) = [];
                            
                            % in the matrix x are all eigenvalues of the matrix which
                            % should be used for MMSE channel estiamtion
                            
                            x = (eigenvalue_t * eigenvalue_h') ./ (eigenvalue_t * eigenvalue_h' + sigma_n2 * ones(length(eigenvalue_t),length(eigenvalue_h)));
                            
                            %here we use rang1 approximation of matrix, it means we are
                            %caulculating eigenvalues of matrix A and B which shoud
                            %approximate MMSE-matrix (kron(B',A) ~ R_h*(R_h + sigma_n2*I)^-1)
                            [U, S, V] = svd(x);
                            b_small = S(1,1)*U(:,1);
                            a_small = V(:,1);
                            
                            b = zeros(14,1);
                            a = zeros(72,1);
                            b(pos_t_large) = b_small;
                            a(pos_h_large) = a_small;
                            
                            A = eigenvector_h*diag(a)*eigenvector_h';
                            B = eigenvector_t*diag(b)*eigenvector_t';
                            
                            H_ls_one_channel_vector = H_ls_one_channel(:);
                            [position_freq,position_time] = find(RefMapping(:,:,tt));
                            if(tt>2)
                                position_freq_matrix = reshape(position_freq,length(position_freq)/2,2);
                                position_time_matrix = reshape(position_time,length(position_time)/2,2);
                                H_ls_one_channel_vector(length(position_time) + 1:end) = [];
                            else
                                position_freq_matrix = reshape(position_freq,length(position_freq)/4,4);
                                position_time_matrix = reshape(position_time,length(position_time)/4,4);
                            end
                            %add "virtual" pilots position, which we will extrapolate
                            if(tt>2)
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1), position_freq_matrix_ext(:,end), position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3, position_time_matrix_ext(:,end) + 6, position_time_matrix_ext(:,end) + 9];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            else
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            end
                            
                            %remove "virtual pilots" which are not nearest to real symbols
                            position_time_ext(position_not_needed) = [];
                            position_freq_ext(position_not_needed) = [];
                            position_to_extrapolate = zeros(length(position_time_ext)-length(position_time),2);
                            extrapolated_values = zeros(length(position_time_ext)-length(position_time),1);
                            three_nearest_points = zeros(length(position_time_ext)-length(position_time),3);
                            i = 1;
                            for pos_i = 1:length(position_time_ext)
                                if(sum(position_freq == position_freq_ext(pos_i) & position_time == position_time_ext(pos_i)))
                                else
                                    % here we find for every "virtual pilot" three
                                    % nearest pilots values and set the value of
                                    % "virtual pilot" as if it would be in plane,
                                    % which is spanned by those nearest points
                                    position_to_extrapolate(i,:) = [position_freq_ext(pos_i),position_time_ext(pos_i)];
                                    [distance , nearest_index] = sort((position_freq - position_freq_ext(pos_i)).^2 + (position_time - position_time_ext(pos_i)).^2);
                                    three_nearest_points(i,:) = nearest_index(1:3);
                                    three_nearest_points_help = sort(three_nearest_points(i,:));
                                    if(three_nearest_points_help == [three_nearest_points_help(2)-1,three_nearest_points_help(2),three_nearest_points_help(2)+1])
                                        positions_help = find(nearest_index < three_nearest_points_help(2) + 6 & nearest_index > three_nearest_points_help(2) - 6);
                                        nearest_index(positions_help) = [];
                                        three_nearest_points(i,3) = nearest_index(1);
                                    end
                                    grad = linsolve([position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,2)) , position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,3)) ; position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,2)) , position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,3))] ,[position_to_extrapolate(i,:)].' - [position_freq(three_nearest_points(i,1)); position_time(three_nearest_points(i,1))]);
                                    extrapolated_values(i) = H_ls_one_channel_vector(three_nearest_points(i,1)) + [H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,2)), H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,3))] * grad;
                                    i = i + 1;
                                end
                            end
                            [XI, YI] = meshgrid(1:LTE_params.Ntot,1:LTE_params.Nsub);
                            channel_estimate_help_ls = griddata([position_freq; position_to_extrapolate(:,1)],[position_time; position_to_extrapolate(:,2)],[H_ls_one_channel_vector; extrapolated_values], XI, YI,channel_interpolation_method);
                            channel_estimate_help_ls = channel_estimate_help_ls.';
                            channel_estimate_help = A*channel_estimate_help_ls*B';
                            channel_estimate(:,:,rr,tt) = channel_estimate_help;
                        end
                        
                        if( strcmp(UE.autocorrelation_matrix_type,'estimated') && UE.realization_num < UE.realization_num_total )
                            %cubic interpolation
                            %add "virtual" pilots position, which we will extrapolate
                            if(tt>2)
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1), position_freq_matrix_ext(:,end), position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3, position_time_matrix_ext(:,end) + 6, position_time_matrix_ext(:,end) + 9];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            else
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            end
                            
                            %remove "virtual pilots" which are not nearest to real symbols
                            position_time_ext(position_not_needed) = [];
                            position_freq_ext(position_not_needed) = [];
                            position_to_extrapolate = zeros(length(position_time_ext)-length(position_time),2);
                            extrapolated_values = zeros(length(position_time_ext)-length(position_time),1);
                            three_nearest_points = zeros(length(position_time_ext)-length(position_time),3);
                            i = 1;
                            for pos_i = 1:length(position_time_ext)
                                if(sum(position_freq == position_freq_ext(pos_i) & position_time == position_time_ext(pos_i)))
                                else
                                    % here we find for every "virtual pilot" three
                                    % nearest pilots values and set the value of
                                    % "virtual pilot" as if it would be in plane,
                                    % which is spanned by those nearest points
                                    position_to_extrapolate(i,:) = [position_freq_ext(pos_i),position_time_ext(pos_i)];
                                    [distance , nearest_index] = sort((position_freq - position_freq_ext(pos_i)).^2 + (position_time - position_time_ext(pos_i)).^2);
                                    three_nearest_points(i,:) = nearest_index(1:3);
                                    three_nearest_points_help = sort(three_nearest_points(i,:));
                                    if(three_nearest_points_help == [three_nearest_points_help(2)-1,three_nearest_points_help(2),three_nearest_points_help(2)+1])
                                        positions_help = find(nearest_index < three_nearest_points_help(2) + 6 & nearest_index > three_nearest_points_help(2) - 6);
                                        nearest_index(positions_help) = [];
                                        three_nearest_points(i,3) = nearest_index(1);
                                    end
                                    grad = linsolve([position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,2)) , position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,3)) ; position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,2)) , position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,3))] ,[position_to_extrapolate(i,:)].' - [position_freq(three_nearest_points(i,1)); position_time(three_nearest_points(i,1))]);
                                    extrapolated_values(i) = H_ls_one_channel_vector(three_nearest_points(i,1)) + [H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,2)), H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,3))] * grad;
                                    i = i + 1;
                                end
                            end
                            [XI, YI] = meshgrid(1:LTE_params.Ntot,1:LTE_params.Nsub);
                            channel_estimate_help = griddata([position_freq; position_to_extrapolate(:,1)],[position_time; position_to_extrapolate(:,2)],[H_ls_one_channel_vector; extrapolated_values], XI, YI,'cubic');
                            channel_estimate(:,:,rr,tt) = channel_estimate_help.';
                            %update of autocorrelation amtrix
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix * UE.realization_num + (channel_estimate(:,:,rr,tt)*channel_estimate(:,:,rr,tt)')/LTE_params.Nsub;
                            UE.realization_num = UE.realization_num + 1;
                            UE.channel_autocorrelation_matrix = UE.channel_autocorrelation_matrix/UE.realization_num;
                            if UE.realization_num == UE.realization_num_total
                                correlation_ceof_vector_help = nan(LTE_params.Ntot,1);
                                for corr_i = 1:LTE_params.Ntot
                                    matrix_help = UE.channel_autocorrelation_matrix(corr_i:LTE_params.Ntot,1:LTE_params.Ntot-corr_i+1);
                                    correlation_ceof_vector_help(corr_i) = mean(diag(matrix_help));
                                end
                                UE.channel_autocorrelation_matrix = toeplitz(correlation_ceof_vector_help);
                            end
                        end
                        interpolation = false;
                        
                    case 'ALMMSE2'
                        H_ls_one_channel = rx_ref_symbols(:,:,tt,rr)./RefSym(:,:,tt);
                        interpolation = false;
                        
                        H_ls_one_channel_vector = H_ls_one_channel(:);
                        [position_freq,position_time] = find(RefMapping(:,:,tt));
                        if(tt>2)
                            position_freq_matrix = reshape(position_freq,length(position_freq)/2,2);
                            position_time_matrix = reshape(position_time,length(position_time)/2,2);
                            H_ls_one_channel_vector(length(position_time) + 1:end) = [];
                        else
                            position_freq_matrix = reshape(position_freq,length(position_freq)/4,4);
                            position_time_matrix = reshape(position_time,length(position_time)/4,4);
                        end
                        %add "virtual" pilots position, which we will extrapolate
                        if(tt>2)
                            position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                            position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                            position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1), position_freq_matrix_ext(:,end), position_freq_matrix_ext(:,end-1) ];
                            position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3, position_time_matrix_ext(:,end) + 6, position_time_matrix_ext(:,end) + 9];
                            position_time_ext = position_time_matrix_ext(:);
                            position_freq_ext = position_freq_matrix_ext(:);
                            position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                        else
                            position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                            position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                            position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1) ];
                            position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3];
                            position_time_ext = position_time_matrix_ext(:);
                            position_freq_ext = position_freq_matrix_ext(:);
                            position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                        end
                        
                        %remove "virtual pilots" which are not nearest to real symbols
                        position_time_ext(position_not_needed) = [];
                        position_freq_ext(position_not_needed) = [];
                        position_to_extrapolate = zeros(length(position_time_ext)-length(position_time),2);
                        extrapolated_values = zeros(length(position_time_ext)-length(position_time),1);
                        three_nearest_points = zeros(length(position_time_ext)-length(position_time),3);
                        i = 1;
                        for pos_i = 1:length(position_time_ext)
                            if(sum(position_freq == position_freq_ext(pos_i) & position_time == position_time_ext(pos_i)))
                            else
                                % here we find for every "virtual pilot" three
                                % nearest pilots values and set the value of
                                % "virtual pilot" as if it would be in plane,
                                % which is spanned by those nearest points
                                position_to_extrapolate(i,:) = [position_freq_ext(pos_i),position_time_ext(pos_i)];
                                [distance , nearest_index] = sort((position_freq - position_freq_ext(pos_i)).^2 + (position_time - position_time_ext(pos_i)).^2);
                                three_nearest_points(i,:) = nearest_index(1:3);
                                three_nearest_points_help = sort(three_nearest_points(i,:));
                                if(three_nearest_points_help == [three_nearest_points_help(2)-1,three_nearest_points_help(2),three_nearest_points_help(2)+1])
                                    positions_help = find(nearest_index < three_nearest_points_help(2) + 6 & nearest_index > three_nearest_points_help(2) - 6);
                                    nearest_index(positions_help) = [];
                                    three_nearest_points(i,3) = nearest_index(1);
                                end
                                grad = linsolve([position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,2)) , position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,3)) ; position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,2)) , position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,3))] ,[position_to_extrapolate(i,:)].' - [position_freq(three_nearest_points(i,1)); position_time(three_nearest_points(i,1))]);
                                extrapolated_values(i) = H_ls_one_channel_vector(three_nearest_points(i,1)) + [H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,2)), H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,3))] * grad;
                                i = i + 1;
                            end
                        end
                        [XI, YI] = meshgrid(1:LTE_params.Ntot,1:LTE_params.Nsub);
                        channel_estimate_help = griddata([position_freq; position_to_extrapolate(:,1)],[position_time; position_to_extrapolate(:,2)],[H_ls_one_channel_vector; extrapolated_values], XI, YI,channel_interpolation_method);
                        
                        time_correlation_temp = channel_estimate_help*channel_estimate_help';
                        freq_correlation_temp = channel_estimate_help.'*conj(channel_estimate_help);
                        
                        nr_of_used_pilots = 2;
                        
                        mapping_tx = RefMapping(:,:,tt);
                        nr_of_pilots = sum(sum(mapping_tx));
                        mapping_tx_double = double(mapping_tx);
                        mapping_tx_double(mapping_tx) = [1:nr_of_pilots];
                        closest_pilots = zeros(LTE_params.Ntot,LTE_params.Nsub,nr_of_used_pilots);
                        
                        max_distance_freq = 0;
                        max_distance_time = 0;
                        
                        for ff = 1:LTE_params.Ntot
                            for ss = 1:LTE_params.Nsub
                                time_distance_temp = (position_time-ss).^2;
                                freq_distance_temp = (position_freq-ff).^2;
                                distance_temp = sqrt(freq_distance_temp+time_distance_temp);
                                [distance_temp_sort index] = sort(distance_temp);
                                closest_pilots(ff,ss,:) = index(1:nr_of_used_pilots);
                                
                                if max_distance_freq < sqrt(freq_distance_temp(index(nr_of_used_pilots)))
                                    max_distance_freq = sqrt(freq_distance_temp(index(nr_of_used_pilots)));
                                end
                                
                                if max_distance_time < sqrt(time_distance_temp(index(nr_of_used_pilots)))
                                    max_distance_time = sqrt(time_distance_temp(index(nr_of_used_pilots)));
                                end
                            end
                        end
                        
                        time_correlation_temp_small = zeros(max_distance_time);
                        freq_correlation_temp_small = zeros(max_distance_freq);
                        
                        time_step = floor(max_distance_time/2);
                        for step_i = 1:ceil((LTE_params.Nsub-max_distance_time)/time_step)
                            time_correlation_temp_small = time_correlation_temp_small + time_correlation_temp((step_i-1)*time_step+1:(step_i-1)*time_step+max_distance_time,(step_i-1)*time_step+1:(step_i-1)*time_step+max_distance_time);
                        end
                        time_correlation_temp_small = time_correlation_temp_small/ceil((LTE_params.Nsub-max_distance_time)/time_step);
                        time_correlation_temp_small = time_correlation_temp_small/LTE_params.Ntot;
                        
                        freq_step = floor(max_distance_freq/2);
                        for step_i = 1:ceil((LTE_params.Ntot-max_distance_freq)/freq_step)
                            freq_correlation_temp_small = freq_correlation_temp_small + freq_correlation_temp((step_i-1)*freq_step+1:(step_i-1)*freq_step+max_distance_freq,(step_i-1)*freq_step+1:(step_i-1)*freq_step+max_distance_freq);
                        end
                        freq_correlation_temp_small = freq_correlation_temp_small/ceil((LTE_params.Ntot-max_distance_freq)/freq_step);
                        freq_correlation_temp_small = freq_correlation_temp_small/LTE_params.Nsub;
                        
                        autocorrelation_matrix_small = kron(time_correlation_temp_small,freq_correlation_temp_small);
                        
                        for ff = 1:LTE_params.Ntot
                            for ss = 1:LTE_params.Nsub
                                
                                
                                
                            end
                        end
                        
                    otherwise
                        error('not supported channel estimation method for fastfading');
                end
                %% channel interpolation
                if(interpolation)
                    switch channel_interpolation_method
                        case {'linear','cubic'}
                            
                            H_ls_one_channel_vector = H_ls_one_channel(:);
                            [position_freq,position_time] = find(RefMapping(:,:,tt));
                            if(tt>2)
                                position_freq_matrix = reshape(position_freq,length(position_freq)/2,2);
                                position_time_matrix = reshape(position_time,length(position_time)/2,2);
                                H_ls_one_channel_vector(length(position_time) + 1:end) = [];
                            else
                                position_freq_matrix = reshape(position_freq,length(position_freq)/4,4);
                                position_time_matrix = reshape(position_time,length(position_time)/4,4);
                            end
                            %add "virtual" pilots position, which we will extrapolate
                            if(tt>2)
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1), position_freq_matrix_ext(:,end), position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3, position_time_matrix_ext(:,end) + 6, position_time_matrix_ext(:,end) + 9];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            else
                                position_freq_matrix_ext = [position_freq_matrix(1,:) - 6;position_freq_matrix; position_freq_matrix(end,:) + 6];
                                position_time_matrix_ext = [position_time_matrix(1,:);position_time_matrix; position_time_matrix(end,:)];
                                position_freq_matrix_ext = [position_freq_matrix_ext(:,2) ,position_freq_matrix_ext, position_freq_matrix_ext(:,end-1) ];
                                position_time_matrix_ext = [position_time_matrix_ext(:,1)-3, position_time_matrix_ext, position_time_matrix_ext(:,end) + 3];
                                position_time_ext = position_time_matrix_ext(:);
                                position_freq_ext = position_freq_matrix_ext(:);
                                position_not_needed = find(position_time_ext < -3 | position_time_ext > LTE_params.Nsub + 3 | position_freq_ext < -6 | position_freq_ext > LTE_params.Ntot + 6);
                            end
                            
                            %remove "virtual pilots" which are not nearest to real symbols
                            position_time_ext(position_not_needed) = [];
                            position_freq_ext(position_not_needed) = [];
                            position_to_extrapolate = zeros(length(position_time_ext)-length(position_time),2);
                            extrapolated_values = zeros(length(position_time_ext)-length(position_time),1);
                            three_nearest_points = zeros(length(position_time_ext)-length(position_time),3);
                            i = 1;
                            for pos_i = 1:length(position_time_ext)
                                if(sum(position_freq == position_freq_ext(pos_i) & position_time == position_time_ext(pos_i)))
                                else
                                    % here we find for every "virtual pilot" three
                                    % nearest pilots values and set the value of
                                    % "virtual pilot" as if it would be in plane,
                                    % which is spanned by those nearest points
                                    position_to_extrapolate(i,:) = [position_freq_ext(pos_i),position_time_ext(pos_i)];
                                    [distance , nearest_index] = sort((position_freq - position_freq_ext(pos_i)).^2 + (position_time - position_time_ext(pos_i)).^2);
                                    three_nearest_points(i,:) = nearest_index(1:3);
                                    three_nearest_points_help = sort(three_nearest_points(i,:));
                                    if(three_nearest_points_help == [three_nearest_points_help(2)-1,three_nearest_points_help(2),three_nearest_points_help(2)+1])
                                        positions_help = find(nearest_index < three_nearest_points_help(2) + 6 & nearest_index > three_nearest_points_help(2) - 6);
                                        nearest_index(positions_help) = [];
                                        three_nearest_points(i,3) = nearest_index(1);
                                    end
                                    grad = linsolve([position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,2)) , position_freq(three_nearest_points(i,1)) - position_freq(three_nearest_points(i,3)) ; position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,2)) , position_time(three_nearest_points(i,1)) - position_time(three_nearest_points(i,3))] ,[position_to_extrapolate(i,:)].' - [position_freq(three_nearest_points(i,1)); position_time(three_nearest_points(i,1))]);
                                    extrapolated_values(i) = H_ls_one_channel_vector(three_nearest_points(i,1)) + [H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,2)), H_ls_one_channel_vector(three_nearest_points(i,1)) - H_ls_one_channel_vector(three_nearest_points(i,3))] * grad;
                                    i = i + 1;
                                end
                            end
                            [XI, YI] = meshgrid(1:LTE_params.Ntot,1:LTE_params.Nsub);
                            channel_estimate_help = griddata([position_freq; position_to_extrapolate(:,1)],[position_time; position_to_extrapolate(:,2)],[H_ls_one_channel_vector; extrapolated_values], XI, YI,channel_interpolation_method);
                            channel_estimate(:,:,rr,tt) = channel_estimate_help.';
                            
                        case 'v4'
                            [position_freq,position_time] = find(RefMapping(:,:,tt));
                            H_ls_one_channel_vector = H_ls_one_channel(:);
                            if(tt>2)
                                H_ls_one_channel_vector(length(position_time) + 1:end) = [];
                            end
                            [XI, YI] = meshgrid(1:LTE_params.Ntot,1:LTE_params.Nsub);
                            channel_estimate_help = griddata(position_freq,position_time,H_ls_one_channel_vector, XI, YI,'v4');
                            channel_estimate(:,:,rr,tt) = channel_estimate_help.';
                        otherwise
                            error('not supported interpolation method for fastfading');
                    end
                end
            otherwise
                error('not supported filtering type')
        end
    end
end

switch ChanMod.filtering
    case 'BlockFading'
        switch channel_estimation_method
            case 'GENIE' %LS for block fading
                %                         A = [];
                %                         for ii = 1:LTE_params.Nsub
                %                             A = [A;eye(ChanMod.nRX)];
                %                         end
                %                         B = [];
                %                         for ii = 1:LTE_params.Nsub
                %                             B = [B,eye(ChanMod.nTX)];
                %                         end
                A = [];
                B = [];
                index_vec = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
                for sub_i = 1:LTE_params.Ntot
                    y_tx_assembled_genie_sub = conj(squeeze(y_tx_assembled_genie(sub_i,index_vec,:)));
                    if ChanMod.nTX == 1
                        y_tx_assembled_genie_sub = y_tx_assembled_genie_sub(:);
                    end
                    if rank(y_tx_assembled_genie_sub)~=ChanMod.nTX
                        A = [A;sub_i];
                    elseif rank(y_tx_assembled_genie_sub)==ChanMod.nTX
                        B = [B;sub_i];
                        y_rx_assembled_sub = conj(squeeze(y_rx_assembled(sub_i,index_vec,:)));
                        y_rx_assembled_sub = y_rx_assembled_sub(:);
                        kron_matrix = kron(eye(ChanMod.nRX),y_tx_assembled_genie_sub);
                        LS_matrix = inv(kron_matrix'*kron_matrix)*kron_matrix';
                        channel_help = LS_matrix*y_rx_assembled_sub;
                        channel_help = reshape(channel_help,ChanMod.nTX,ChanMod.nRX)';
                        channel_estimate(sub_i,1,:,:) = channel_help;
                    end
                    %y_tx_assembled_genie_sub =
                    %y_tx_assembled_genie_sub(:);
                end
                
                for tt = 1:ChanMod.nTX
                    for rr = 1:ChanMod.nRX
                        channel_estimate_help = interp1(B,squeeze(channel_estimate(B,1,rr,tt)),1:LTE_params.Ntot,channel_interpolation_method,'extrap');
                        channel_estimate_help = channel_estimate_help.';
                        channel_estimate(:,1,rr,tt) = channel_estimate_help;
                    end
                end
                channel_estimate = repmat(channel_estimate(:,1,:,:),[1,LTE_params.Nsub,1,1]);
        end
end