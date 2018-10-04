function H_est_large = LTE_ICI_estimator(LTE_params,H_est,ChanMod,ChanMod_output)
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
% date of creation: 2011/01/27
% last changes: 2011/01/27  Simko

%% ICI estimation
    H_est = H_est + sqrt(LTE_params.channel_noise)*(randn(size(H_est))+sqrt(-1)*randn(size(H_est)))/sqrt(2);
    H_est_large = zeros(LTE_params.Ntot,LTE_params.Ntot,LTE_params.Nsub,ChanMod.nRX,ChanMod.nTX);
    for tt = 1:ChanMod.nTX
        for rr = 1:ChanMod.nRX
            switch LTE_params.ICI_est_type
                case 'PERFECT'
                    H_est_large(:,:,:,rr,tt) = ChanMod_output.genie.H_fft_matrix(:,:,:,rr,tt);
                case '2nd'
                    %% higher order method 2nd type
                    order_of_approximation = LTE_params.order;
                    channel_samples_position = [LTE_params.Ntot/2];
                    for ss = 2:LTE_params.Nsub
                        channel_samples_position = [channel_samples_position;channel_samples_position(end) + LTE_params.Ntot];
                    end
                    
                    
                    
                    
                    for ss = 1:LTE_params.Nsub
                        
                        channel_samples_position_used = channel_samples_position - channel_samples_position(ss);
                        A = [ones(LTE_params.Nsub,1)];
                        
                        A_long = [ones(LTE_params.Ntot,1)];
                        
                        for o_i = 1:order_of_approximation
                            A = [A,channel_samples_position_used.^o_i];
                            
                            start = -LTE_params.Ntot/2;
                            stop = start + LTE_params.Ntot - 1;
                            A_long_vec = [start:stop]';
                            A_long = [A_long,A_long_vec.^o_i];
                        end
                        
                        
                        %gram schmidt
                        [m n] = size(A);
                        Q = zeros(m,n);
                        R = zeros(m,n);
                        
                        Q_long = zeros(size(A_long));
                        
                        for j=1:n
                            v = A(:,j);
                            
                            v_long = A_long(:,j);
                            for i=1:j-1
                                R(i,j) = Q(:,i)'*A(:,j);
                                v = v - R(i,j)*Q(:,i);
                                
                                v_long = v_long - R(i,j)*Q_long(:,i);
                            end
                            R(j,j)=norm(v);
                            Q(:,j)=v/R(j,j);
                            
                            Q_long(:,j)=v_long/R(j,j);
                        end
                    
                        LS_matrix = Q';
                        H_temp = squeeze(H_est(:,:,rr,tt));
                        channel_coef = LS_matrix*H_temp.';
                        
                        
                        
                        for oo = 0:order_of_approximation
                                                       
                            DFT_matrix = sqrt(1/LTE_params.Ntot)*dftmtx(LTE_params.Ntot);
                            
                            ici_temp_small = diag(channel_coef(oo+1,:))*DFT_matrix*diag(Q_long(:,oo+1))*DFT_matrix';
                            H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                            
                        end
                        ici_est_temp = squeeze(H_est_large(:,:,ss,rr,tt));
                        ici_est_temp(logical(eye(72))) = H_est(:,ss,rr,tt);
                        H_est_large(:,:,ss,rr,tt) = ici_est_temp;
                    end
                case 'thomas'
                    %% higher order method 2nd type
                    order_of_approximation = LTE_params.order;
                    channel_samples_position = [LTE_params.Ntot/2];
                    for ss = 2:LTE_params.Nsub
                        channel_samples_position = [channel_samples_position;channel_samples_position(end) + LTE_params.Ntot];
                    end
                    
                    c = LTE_params.speed_of_light;
                    f = LTE_params.carrier_freq;  % Frequency at which our system operates
                    v = UE.user_speed;  %speed at which we move
                    f_d = v*f/c;   % maximum radian Doppler frequency
                    nuM = f_d*LTE_params.Tsubframe;
                    
                    %*** set up inital basisfunctions for channel interpolation
                    [Ubi,Sbi] = dpss(LTE_params.Nsub*LTE_params.Ntot,nuM,LTE_params.order+1);
                    
                    
                    %"pilot" position
                    P = channel_samples_position;
                    A = Ubi(P,:);
                    
                    for ss = 1:LTE_params.Nsub
                        
                        LS_matrix = pinv(A);
                        %LS_matrix = Q';
                        H_temp = squeeze(H_est(:,:,rr,tt));
                        channel_coef = LS_matrix*H_temp.';
                        
                        start = LTE_params.Ntot*(ss-1)+1;
                        stop =  LTE_params.Ntot*ss;
                        for oo = 1:order_of_approximation
                            
                            time_matrix = diag(Ubi(start:stop,oo+1));
                            
                            DFT_matrix = sqrt(1/LTE_params.Ntot)*dftmtx(LTE_params.Ntot);
                            
                            ici_temp_small = diag(channel_coef(oo+1,:))*DFT_matrix*(time_matrix)*DFT_matrix';
                            H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                            
                        end
%                         ici_est_temp = squeeze(H_est_large(:,:,ss,rr,tt));
%                         ici_est_temp(logical(eye(72))) = H_est(:,ss,rr,tt);
%                         H_est_large(:,:,ss,rr,tt) = ici_est_temp;
                    end
                case 'thomas_2nd'
                    %% higher order method 2nd type
                    order_of_approximation = LTE_params.order;
                    channel_samples_position = [LTE_params.Ntot/2];
                    for ss = 2:LTE_params.Nsub
                        channel_samples_position = [channel_samples_position;channel_samples_position(end) + LTE_params.Ntot];
                    end
                    
                    c = LTE_params.speed_of_light;
                    f = LTE_params.carrier_freq;  % Frequency at which our system operates
                    v = UE.user_speed;  %speed at which we move
                    f_d = v*f/c;   % maximum radian Doppler frequency
                    nuM = f_d*LTE_params.Tsubframe;
                    
                    %*** set up inital basisfunctions for channel interpolation
                    [Ubi,Sbi] = dpss(LTE_params.Nsub*LTE_params.Ntot,nuM,LTE_params.order+1);
                    
                    
                    %"pilot" position
                    P = channel_samples_position;
                    A = Ubi(P,:);
                    
                    
                    %gram schmidt
                    [m n] = size(A);
                    Q = zeros(m,n);
                    R = zeros(m,n);
                    
                    Q_long = zeros(size(Ubi));
                    
                    for j=1:n
                        v = A(:,j);
                        
                        v_long = Ubi(:,j);
                        for i=1:j-1
                            R(i,j) = Q(:,i)'*A(:,j);
                            v = v - R(i,j)*Q(:,i);
                            
                            v_long = v_long - R(i,j)*Q_long(:,i);
                        end
                        R(j,j)=norm(v);
                        Q(:,j)=v/R(j,j);
                        
                        Q_long(:,j)=v_long/R(j,j);
                    end
                    
                    for ss = 1:LTE_params.Nsub
                        
                        %LS_matrix = pinv(A);
                        LS_matrix = Q';
                        H_temp = squeeze(H_est(:,:,rr,tt));
                        channel_coef = LS_matrix*H_temp.';
                        
                        start = LTE_params.Ntot*(ss-1)+1;
                        stop =  LTE_params.Ntot*ss;
                        for oo = 1:order_of_approximation
                            
                            time_matrix = diag(Q_long(start:stop,oo+1));
                            
                            DFT_matrix = sqrt(1/LTE_params.Ntot)*dftmtx(LTE_params.Ntot);
                            
                            ici_temp_small = diag(channel_coef(oo+1,:))*DFT_matrix*(time_matrix)*DFT_matrix';
                            H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                            
                        end
                        ici_est_temp = squeeze(H_est_large(:,:,ss,rr,tt));
                        ici_est_temp(logical(eye(72))) = H_est(:,ss,rr,tt);
                        H_est_large(:,:,ss,rr,tt) = ici_est_temp;
                    end
                case 'thomas_1st'
                    %% higher order method
                    channel_samples_position = [LTE_params.NfftCP{1}/2];
                    order_of_approximation = LTE_params.order;
                    table_start_end = nan(LTE_params.Nsub,2);
                    table_start_end(1,1) = 1;
                    table_start_end(1,2) = LTE_params.NfftCP{1};
                    
                    for ss = 2:LTE_params.Nsub
                        if ss<7
                            channel_samples_position = [channel_samples_position;LTE_params.NfftCP{1} + (ss-2)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{2}/2)];  % here its assumed, that the samples are in the middle of an ofdm symbol, but it can be easily changed
                            table_start_end(ss,1) = table_start_end(ss-1,2) + 1;
                            table_start_end(ss,2) = table_start_end(ss,1) + LTE_params.NfftCP{2}-1;
                        elseif ss==7
                            channel_samples_position = [channel_samples_position;LTE_params.NfftCP{1} + (ss-2)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{1}/2)];  % here its assumed, that the samples are in the middle of an ofdm symbol, but it can be easily changed
                            table_start_end(ss,1) = table_start_end(ss-1,2) + 1;
                            table_start_end(ss,2) = table_start_end(ss,1) + LTE_params.NfftCP{1}-1;
                        else
                            channel_samples_position = [channel_samples_position;2*LTE_params.NfftCP{1} + (ss-3)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{2}/2)];
                            table_start_end(ss,1) = table_start_end(ss-1,2) + 1;
                            table_start_end(ss,2) = table_start_end(ss,1) + LTE_params.NfftCP{2}-1;
                        end
                    end
                    
                    c = LTE_params.speed_of_light;
                    f = LTE_params.carrier_freq;  % Frequency at which our system operates
                    v = LTE_params.UE_config.user_speed;  %speed at which we move
                    f_d = v*f/c;   % maximum radian Doppler frequency
                    %nuM = f_d*LTE_params.Tsubframe;
                    nM=LTE_params.NfftCP{1}*2+LTE_params.NfftCP{2}*12;
                    nuM = f_d*LTE_params.SamplingTime*nM;
                    %nuM = f_d*LTE_params.SamplingTime*137;
                    
                    %*** set up inital basisfunctions for channel interpolation
                    [Ubi,Sbi] = dpss(nM,nuM,LTE_params.order+1);
                    
                    
                    %"pilot" position
                    P = channel_samples_position;
                    A = Ubi(P,:);                    
                                                                  
                    LS_matrix = pinv(A);
                    H_temp = squeeze(H_est(:,:,rr,tt));
                    channel_coef = LS_matrix*H_temp.';
                        
                    for ss = 1:LTE_params.Nsub


                        
                        for oo = 0:order_of_approximation
                            
                            
                            
                            if(ss == 1 || ss == 7) %the calculation of the output for longer symbols
                                start = table_start_end(ss,1);
                                stop =table_start_end(ss,2);
                                time_matrix = diag(Ubi(start:stop,oo+1));
                                
                                CP_length = LTE_params.NfftCP{1} - LTE_params.Nfft;
                                F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                
                                ici_temp_small = DFT_matrix*F_CP_rem*(time_matrix)*F_CP_add*DFT_matrix';
                                ici_temp_small = diag(channel_coef(oo+1,:))*ici_temp_small([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                                
                            else %the calculation of the output for shorter symbols
                                start = table_start_end(ss,1);
                                stop =table_start_end(ss,2);
                                time_matrix = diag(Ubi(start:stop,oo+1));
                                
                                CP_length = LTE_params.NfftCP{2} - LTE_params.Nfft;
                                F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                
                                ici_temp_small = DFT_matrix*F_CP_rem*(time_matrix)*F_CP_add*DFT_matrix';
                                ici_temp_small = diag(channel_coef(oo+1,:))*ici_temp_small([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                                
                            end
                        end
                           %ici_est_temp = squeeze(H_est_large(:,:,ss,rr,tt));
                          %ici_est_temp(logical(eye(72))) = H_est(:,ss,rr,tt);
                           %H_est_large(:,:,ss,rr,tt) = ici_est_temp;
                    end
                case '1st'
                    %% higher order method
                    channel_samples_position = [LTE_params.NfftCP{1}/2];
                    order_of_approximation = LTE_params.order;
                    for ss = 2:LTE_params.Nsub
                        if ss<7
                            channel_samples_position = [channel_samples_position;LTE_params.NfftCP{1} + (ss-2)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{2}/2)];  % here its assumed, that the samples are in the middle of an ofdm symbol, but it can be easily changed
                        elseif ss==7
                            channel_samples_position = [channel_samples_position;LTE_params.NfftCP{1} + (ss-2)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{1}/2)];  % here its assumed, that the samples are in the middle of an ofdm symbol, but it can be easily changed
                        else
                            channel_samples_position = [channel_samples_position;2*LTE_params.NfftCP{1} + (ss-3)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{2}/2)];
                        end
                    end
                    
                    
                    
                    
                    for ss = 1:LTE_params.Nsub
                        
                        channel_samples_position_used = channel_samples_position - channel_samples_position(ss);
                        A = [ones(LTE_params.Nsub,1)];
                        for o_i = 1:order_of_approximation
                            A = [A,channel_samples_position_used.^o_i];
                        end
                        
                        
                        LS_matrix = pinv(A);
                        H_temp = squeeze(H_est(:,:,rr,tt));
                        channel_coef = LS_matrix*H_temp.';
                        
                        for oo = 0:order_of_approximation
                            
                            
                            
                            
                            
                            if(ss == 1 || ss == 7) %the calculation of the output for longer symbols
                                start = -LTE_params.NfftCP{1}/2;
                                stop = start + LTE_params.NfftCP{1} - 1;
                                time_matrix = diag(start:stop);
                                
                                CP_length = LTE_params.NfftCP{1} - LTE_params.Nfft;
                                F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                
                                ici_temp_small = DFT_matrix*F_CP_rem*(time_matrix.^oo)*F_CP_add*DFT_matrix';
                                ici_temp_small = diag(channel_coef(oo+1,:))*ici_temp_small([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                                
                            else %the calculation of the output for shorter symbols
                                start = -ceil(LTE_params.NfftCP{1}/2);
                                stop = start + LTE_params.NfftCP{2} - 1;
                                time_matrix = diag(start:stop);
                                
                                CP_length = LTE_params.NfftCP{2} - LTE_params.Nfft;
                                F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                
                                ici_temp_small = DFT_matrix*F_CP_rem*(time_matrix.^oo)*F_CP_add*DFT_matrix';
                                ici_temp_small = diag(channel_coef(oo+1,:))*ici_temp_small([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                                
                            end
                        end
                        ici_est_temp = squeeze(H_est_large(:,:,ss,rr,tt));
                        ici_est_temp(logical(eye(72))) = H_est(:,ss,rr,tt);
                        H_est_large(:,:,ss,rr,tt) = ici_est_temp;
                    end
                case '1st_ort'
                    %% higher order method
                    channel_samples_position = [LTE_params.NfftCP{1}/2];
                    order_of_approximation = LTE_params.order;
                    table_start_end = nan(LTE_params.Nsub,2);
                    table_start_end(1,1) = 1;
                    table_start_end(1,2) = LTE_params.NfftCP{1};
                    for ss = 2:LTE_params.Nsub
                        if ss<7
                            channel_samples_position = [channel_samples_position;LTE_params.NfftCP{1} + (ss-2)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{2}/2)];  % here its assumed, that the samples are in the middle of an ofdm symbol, but it can be easily changed
                            table_start_end(ss,1) = table_start_end(ss-1,2) + 1;
                            table_start_end(ss,2) = table_start_end(ss,1) + LTE_params.NfftCP{2}-1;
                        elseif ss==7
                            channel_samples_position = [channel_samples_position;LTE_params.NfftCP{1} + (ss-2)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{1}/2)];  % here its assumed, that the samples are in the middle of an ofdm symbol, but it can be easily changed
                            table_start_end(ss,1) = table_start_end(ss-1,2) + 1;
                            table_start_end(ss,2) = table_start_end(ss,1) + LTE_params.NfftCP{1}-1;
                        else
                            channel_samples_position = [channel_samples_position;2*LTE_params.NfftCP{1} + (ss-3)*LTE_params.NfftCP{2} + ceil(LTE_params.NfftCP{2}/2)];
                            table_start_end(ss,1) = table_start_end(ss-1,2) + 1;
                            table_start_end(ss,2) = table_start_end(ss,1) + LTE_params.NfftCP{2}-1;
                        end
                    end
                    
                    %channel_samples_position = ChanMod_output.genie.channel_samples_position(:,rr,tt);
                    
                    A_long = [ones(LTE_params.TxSymbols,1)];
                    A = [ones(LTE_params.Nsub,1)];
                    
                    for o_i = 1:order_of_approximation
                        A = [A,channel_samples_position.^o_i];
                        
                        A_long_vec = [1:LTE_params.TxSymbols]';
                        A_long = [A_long,A_long_vec.^o_i];
                    end
                    
                    
                    %gram schmidt
                    [m n] = size(A);
                    Q = zeros(m,n);
                    R = zeros(m,n);
                    
                    Q_long = zeros(size(A_long));
                    
                    for j=1:n
                        v = A(:,j);
                        
                        v_long = A_long(:,j);
                        for i=1:j-1
                            R(i,j) = Q(:,i)'*A(:,j);
                            v = v - R(i,j)*Q(:,i);
                            
                            v_long = v_long - R(i,j)*Q_long(:,i);
                        end
                        R(j,j)=norm(v);
                        Q(:,j)=v/R(j,j);
                        
                        Q_long(:,j)=v_long/R(j,j);
                    end
                    
                    LS_matrix = Q';
                    H_temp = squeeze(H_est(:,:,rr,tt));
                    channel_coef = LS_matrix*H_temp.';
                    
                    for ss = 1:LTE_params.Nsub
                                              
                        for oo = 0:order_of_approximation
                                                        
                            if(ss == 1 || ss == 7) %the calculation of the output for longer symbols
                                start = table_start_end(ss,1);
                                stop =table_start_end(ss,2);
                                time_matrix = diag(Q_long(start:stop,oo+1));
                                
                                CP_length = LTE_params.NfftCP{1} - LTE_params.Nfft;
                                F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                
                                ici_temp_small = DFT_matrix*F_CP_rem*time_matrix*F_CP_add*DFT_matrix';
                                ici_temp_small = diag(channel_coef(oo+1,:))*ici_temp_small([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                                
                            else %the calculation of the output for shorter symbols
                                start = table_start_end(ss,1);
                                stop = table_start_end(ss,2);
                                time_matrix = diag(Q_long(start:stop,oo+1));
                                
                                CP_length = LTE_params.NfftCP{2} - LTE_params.Nfft;
                                F_CP_add = [[zeros(CP_length,LTE_params.Nfft-CP_length), eye(CP_length)];eye(LTE_params.Nfft)];
                                F_CP_rem = [zeros(LTE_params.Nfft,CP_length),eye(LTE_params.Nfft)];
                                DFT_matrix = sqrt(1/LTE_params.Nfft)*dftmtx(LTE_params.Nfft);
                                
                                ici_temp_small = DFT_matrix*F_CP_rem*time_matrix*F_CP_add*DFT_matrix';
                                ici_temp_small = diag(channel_coef(oo+1,:))*ici_temp_small([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1],[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]);
                                H_est_large(:,:,ss,rr,tt) = H_est_large(:,:,ss,rr,tt) + ici_temp_small;
                                
                            end
                        end
                        %ici_est_temp = squeeze(H_est_large(:,:,ss,rr,tt));
                        %ici_est_temp(logical(eye(72))) = H_est(:,ss,rr,tt);
                        %H_est_large(:,:,ss,rr,tt) = ici_est_temp;
                    end
                case 'linear'
                    %% linear method
                    for ss = 1:LTE_params.Nsub
                        one_before = ss - 1;
                        if one_before == 0
                            one_before = nan;
                        end
                        one_after = ss + 1;
                        if one_after == LTE_params.Nsub+1
                            one_after = nan;
                        end
                        
                        if ss == 1
                            delta_next = squeeze(H_est(:,one_after,rr,tt)) - squeeze(H_est(:,ss,rr,tt));
                            delta_next = delta_next/LTE_params.Ntot;
                            DFT_matrix = sqrt(1/LTE_params.Ntot)*dftmtx(LTE_params.Ntot);
                            
                            mean_point = ceil(LTE_params.Ntot/2);
                            vec_next = zeros(LTE_params.Ntot,1);
                            vec_next = -(mean_point-1):LTE_params.Ntot-mean_point;
                            mat_next = diag(vec_next);
                            ici_est_temp = diag(delta_next) * DFT_matrix * mat_next * DFT_matrix';
                        elseif ss == LTE_params.Nsub
                            delta_previous = squeeze(H_est(:,one_before,rr,tt)) - squeeze(H_est(:,ss,rr,tt));
                            delta_previous = delta_previous/LTE_params.Ntot;
                            DFT_matrix = sqrt(1/LTE_params.Ntot)*dftmtx(LTE_params.Ntot);
                            
                            mean_point = ceil(LTE_params.Ntot/2);
                            vec_previous = zeros(LTE_params.Ntot,1);
                            vec_previous = -(mean_point-1):LTE_params.Ntot-mean_point;
                            mat_previous = diag(vec_previous);
                            ici_est_temp = diag(delta_previous) * DFT_matrix * mat_previous * DFT_matrix';
                        else
                            delta_previous = squeeze(H_est(:,ss,rr,tt)) - squeeze(H_est(:,one_before,rr,tt));
                            delta_previous = delta_previous/LTE_params.Ntot;
                            delta_next = squeeze(H_est(:,one_after,rr,tt)) - squeeze(H_est(:,ss,rr,tt));
                            delta_next = delta_next/LTE_params.Ntot;
                            DFT_matrix = sqrt(1/LTE_params.Ntot)*dftmtx(LTE_params.Ntot);
                            
                            mean_point = ceil(LTE_params.Ntot/2);
                            vec_previous = zeros(LTE_params.Ntot,1);
                            vec_next = zeros(LTE_params.Ntot,1);
                            vec_previous(1:mean_point) = -(mean_point-1):0;
                            vec_next(mean_point:end) = 0:LTE_params.Ntot-mean_point;
                            mat_previous = diag(vec_previous);
                            mat_next = diag(vec_next);
                            ici_est_temp = diag(delta_previous) * DFT_matrix * mat_previous * DFT_matrix' + diag(delta_next) * DFT_matrix * mat_next * DFT_matrix';
                        end
                        
                        ici_est_temp(logical(eye(72))) = H_est(:,ss,rr,tt);
                        H_est_large(:,:,ss,rr,tt) = ici_est_temp;
                        
                    end
            end
        end
    end