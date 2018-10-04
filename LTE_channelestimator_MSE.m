function MSE_tot = LTE_channelestimator_MSE(sigma_n2,autocorrelation_tmp,LTE_params,nTX,speed)

RefMapping_tot = LTE_params.Reference_Signal(1,1).RefMapping;

switch  LTE_params.UE_config.channel_estimation_method
    case 'MMSE'
        switch LTE_params.ChanMod_config.filtering
            case 'BlockFading'
                MSE_tot = zeros(nTX,LTE_params.Ntot);
                for ref_count = 1:size(RefMapping_tot,3)
                    switch ref_count
                        case {1,2}
                            noise_reduce = 2;
                        case {3,4}
                            noise_reduce = 1;
                    end
                    ref_numb = LTE_params.Ntot/3;
                    if ref_count == 1
                        RefMapping = RefMapping_tot(:,:,ref_count);
                        [position_freq,position_time] = find(RefMapping(:,:));
                        position = reshape(position_freq,length(position_freq)/2,2);
                        permutation_matrix = zeros(LTE_params.Ntot);
                        permutation_vector = 1:LTE_params.Ntot;
                        permutation_vector(position(:,1)) = [];
                        permutation_vector = [sort(position(:,1));permutation_vector'];
                        for pos_i = 1:LTE_params.Ntot
                            permutation_matrix(pos_i,permutation_vector(pos_i)) = 1;
                        end
                        perm_tmp = zeros(size(permutation_matrix));
                        perm_tmp2 = zeros(size(permutation_matrix));
                        perm_tmp(ref_numb+1:end,:) = permutation_matrix(ref_numb+1:end,:);
                        [C1,pos_temp] = find(perm_tmp);
                        perm_tmp2(1:ref_numb,:) = permutation_matrix(1:ref_numb,:);
                        [C1,pos_temp2] = find(perm_tmp2);
                        pos_temp = [pos_temp2;pos_temp];
                    end
                    noise = sigma_n2*eye(ref_numb);
                    autocorrelation = permutation_matrix * autocorrelation_tmp * permutation_matrix';
                    R_hLS = autocorrelation(:,1:ref_numb);
                    R_LS = autocorrelation(1:ref_numb,1:ref_numb);
                    A = R_hLS/(R_LS+noise);
                    MSE = zeros(LTE_params.Ntot,1);
                    for ind1 = 1:length(pos_temp)
                        C= (A(ind1,:)'*A(ind1,:));
                        MSE(ind1,:) = 1-2*real(A(ind1,:)*autocorrelation(ind1,1:ref_numb)')+abs(sum(sum(C.*((autocorrelation(1:ref_numb,1:ref_numb).')+sigma_n2/noise_reduce*eye(ref_numb)))));
                    end
                    MSE = permutation_matrix'*MSE;
                    MSE_tot(ref_count,:) = MSE;
                end
%                 plot(MSE_tot(1,:),'g')
%                 hold on
%                 grid on
            case 'FastFading'
                MSE_tot = zeros(nTX,LTE_params.Ntot,LTE_params.Nsub);
                for tt = 1:nTX
                    [position_freq,position_time] = find(RefMapping_tot(:,:,tt));
                    if(tt>2)
                        ref_numb = length(position_freq);
                        position_freq_matrix = reshape(position_freq,ref_numb/2,2);
                        position_time_matrix = reshape(position_time,ref_numb/2,2);
                    else
                        ref_numb = length(position_freq);
                        position_freq_matrix = reshape(position_freq,ref_numb/4,4);
                        position_time_matrix = reshape(position_time,ref_numb/4,4);
                    end
                    c = LTE_params.speed_of_light;
                    f = LTE_params.carrier_freq;  % Frequency at which our system operates
                    v = speed;  %speed at which we move
                    w_d = 2*pi*v*f/c;   % maximum radian Doppler frequency
                    time_autocorrelation = toeplitz(besselj(0,w_d*(0:13)*137*LTE_params.SamplingTime));
                    channel_autocorrelation_full = kron(time_autocorrelation,autocorrelation_tmp);
                    permutation_matrix = zeros(LTE_params.Ntot*LTE_params.Nsub);
                    permutation_vector = 1:(LTE_params.Ntot*LTE_params.Nsub);
                    position = position_freq_matrix + LTE_params.Ntot*(position_time_matrix-1);
                    position = position(:);
                    permutation_vector(position) = [];
                    permutation_vector = [position;permutation_vector'];
                    for pos_i = 1:(LTE_params.Ntot*LTE_params.Nsub)
                        permutation_matrix(pos_i,permutation_vector(pos_i)) = 1;
                    end
                    perm_tmp = zeros(size(permutation_matrix));
                    perm_tmp2 = zeros(size(permutation_matrix));
                    perm_tmp(ref_numb+1:end,:) = permutation_matrix(ref_numb+1:end,:);
                    [C1,pos_temp] = find(perm_tmp);
                    perm_tmp2(1:ref_numb,:) = permutation_matrix(1:ref_numb,:);
                    [C1,pos_temp2] = find(perm_tmp2);
                    pos_temp = [pos_temp2;pos_temp];
                    noise = sigma_n2*eye(ref_numb);
                    autocorrelation = permutation_matrix * channel_autocorrelation_full * permutation_matrix';
                    R_hLS = autocorrelation(:,1:ref_numb);
                    R_LS = autocorrelation(1:ref_numb,1:ref_numb);
                    A = R_hLS/(R_LS+noise);
                    MSE = zeros(LTE_params.Ntot*LTE_params.Nsub,1);
                    for ind1 = 1:length(pos_temp)
                        C= (A(ind1,:)'*A(ind1,:));
                        MSE(ind1,:) = 1-2*real(A(ind1,:)*autocorrelation(ind1,1:ref_numb)')+abs(sum(sum(C.*((autocorrelation(1:ref_numb,1:ref_numb).')+sigma_n2*eye(ref_numb)))));
                    end
                    MSE_temp = permutation_matrix'*MSE;
                    MSE_tot(tt,:,:) = reshape(MSE_temp,LTE_params.Ntot,LTE_params.Nsub);
                end
%                 surf(squeeze(MSE_tot(1,:,:)))  
%                 figure(2)
%                 image(RefMapping_tot(:,:,tt)*100)
        end
                
    case 'LS'
        MSE_tot = zeros(nTX,LTE_params.Ntot);
        for ref_count = 1:size(RefMapping_tot,3)
            switch ref_count
                case {1,2}
                    noise_reduce = 2;
                    RefMapping = logical(RefMapping_tot(:,1,ref_count)) | logical(RefMapping_tot(:,5,ref_count));
                case {3,4}
                    noise_reduce = 1;
                    RefMapping = logical(RefMapping_tot(:,2,ref_count)) | logical(RefMapping_tot(:,9,ref_count));
            end
            ref_numb = LTE_params.Ntot/3;
            [position_freq,position_time] = find(RefMapping,1,'first');
            [position_freq2,position_time2] = find(RefMapping,1,'last');
            deltakf = 3;
            plow = position_freq;
            phigh = position_freq+deltakf;
            MSE_tmp = zeros(1,LTE_params.Ntot);
            for d1 = 1:LTE_params.Ntot
                if d1 < position_freq      %% Extrapolation
                    c1 = -1/deltakf;
                    c2 = 1-c1;
                    MSE_tmp(d1) = c1^2+c2^2+1+(c1^2+c2^2)/noise_reduce*sigma_n2+2*c1*c2*real(autocorrelation_tmp(plow,phigh))-2*c2*real(autocorrelation_tmp(d1,plow))-2*c1*real(autocorrelation_tmp(d1,phigh));                    
                elseif d1 > position_freq2
                    c1 = 1/deltakf;
                    c2 = -c1;
                    c1 = 1+c1;
                    MSE_tmp(d1) = c1^2+c2^2+1+(c1^2+c2^2)/noise_reduce*sigma_n2+2*c1*c2*real(autocorrelation_tmp(plow,phigh))-2*c2*real(autocorrelation_tmp(d1,plow))-2*c1*real(autocorrelation_tmp(d1,phigh));
                else
                    if d1 == phigh+1
                        plow = plow+deltakf;
                        phigh = phigh+deltakf;
                    end 
                    c1 = (d1-plow)/deltakf;
                    c2 = 1-c1;      
                    MSE_tmp(d1) = c1^2+c2^2+1+(c1^2+c2^2)/noise_reduce*sigma_n2+2*c1*c2*real(autocorrelation_tmp(plow,phigh))-2*c2*real(autocorrelation_tmp(d1,plow))-2*c1*real(autocorrelation_tmp(d1,phigh));
                end
            end
            MSE_tot(ref_count,:) = MSE_tmp;    
%             for deltad1 = 0:3
%                 c1 = deltad1/deltakf;
%                 c2 = 1-c1;
%                 MSE_tmp(deltad1+1) = c1^2+c2^2+1+(c1^2+c2^2)/noise_reduce*sigma_n2+2*c1*c2*real(autocorrelation_tmp(1,deltakf))-2*c1*real(autocorrelation_tmp(1,deltakf))-2*c2*sum(PDP_lin.*cos(2*pi*fs*deltad1*PDP_dB(2,:)));
%             end
        end
    case 'PERFECT'
        MSE_tot = zeros(nTX,LTE_params.Ntot);
    otherwise
        MSE_tot = zeros(nTX,LTE_params.Ntot);
end
 
% mean(MSE_tot,2)
% plot(MSE_tot(2,:),'r')
% plot(MSE_tot(3,:),'b')
% % plot(MSE_tot(4,:),'k')