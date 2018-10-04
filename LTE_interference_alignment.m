% Interference Alignment Simulator
% Various IA algorithms
% Author: Joerg Reitterer, jreitter@nt.tuwien.ac.at
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

function [V, U, channel_error] = LTE_interference_alignment(ia_mode, d, P, IA_thresh, IA_max_iterations, LTE_params, ChanMod_output, BS, UE, sigma_n2)

% LTE_params.nBS = 3;
K = LTE_params.nBS;
% nRX = UE{1,1}.nRX;
% nTX = BS{1,1}.nTX;
nRX = LTE_params.UE_config.nRX;
nTX = LTE_params.BS_config.nTx; % why is the "x" in "nTx" not uppercase? typo?

V = cell(K,LTE_params.Ntot,LTE_params.Nsub);
U = cell(K,LTE_params.Ntot,LTE_params.Nsub);

% Feasibility condition for symmetric systems according to
% C. M. Yetis, S. A. Jafar and A. H. Kayran,
% "Feasibility conditions for interference alignment," 
% in IEEE Global Telecommunications Conference 2009, GLOBECOM'09, pp. 1-6, November 2009.
if (nTX*nRX-(K+1)*d(1) < 0)
    warning('System not proper, i.e. almost surely not feasible!');
end

H_fft_noisy = cell(LTE_params.nBS);
for i = 1:K
    for j = 1:K
        H_fft_noisy{i,j} = ChanMod_output{i,j}.genie.H_fft + sqrt(mean(mean(mean(mean(abs(ChanMod_output{i,j}.genie.H_fft).^2))))) * 10^(-LTE_params.IA_sigma_H2_E2_ratio/20) * 1/sqrt(2)*(randn(size(ChanMod_output{i,j}.genie.H_fft)) + 1i*randn(size(ChanMod_output{i,j}.genie.H_fft)));
    end
end

switch ia_mode
    case 'closed_form'
        if (K ~= 3)
            error('Closed form IA solution does not exist for K ~= 3.');
        end

        if (nRX ~= nTX)
            error('Closed form IA solution does not exist for nRX ~= nTX.');
        end
            for m = 1:LTE_params.Nsub/LTE_params.IA_time_granularity
                time_index = LTE_params.IA_time_granularity*(m-1)+1:LTE_params.IA_time_granularity*m;
                for n = 1:LTE_params.Ntot/LTE_params.IA_freq_granularity
                    freq_index = LTE_params.IA_freq_granularity*(n-1)+1:LTE_params.IA_freq_granularity*n;
                    H = get_H(H_fft_noisy,freq_index,time_index(1),LTE_params);
                    E = (H{3,1} \ H{3,2}) * (H{1,2} \ H{1,3}) * (H{2,3} \ H{2,1});

                    [A,~] = eigs(E);

                    %% Compute Precoding Matrices
                    V(1,freq_index,time_index) = repcell(A(:,1:d(1)),[1 LTE_params.IA_freq_granularity LTE_params.IA_time_granularity]);
                    V(2,freq_index,time_index) = repcell((H{3,2} \ H{3,1}) * V{1,freq_index(1),time_index(1)},[1 LTE_params.IA_freq_granularity LTE_params.IA_time_granularity]);
                    V(3,freq_index,time_index) = repcell((H{2,3} \ H{2,1}) * V{1,freq_index(1),time_index(1)},[1 LTE_params.IA_freq_granularity LTE_params.IA_time_granularity]);

                    %% Compute Interference Suppression Matrices
                    U(1,freq_index,time_index) = repcell(null((H{1,2}*V{2,freq_index(1),time_index(1)})'),[1 LTE_params.IA_freq_granularity LTE_params.IA_time_granularity]);
                    U(2,freq_index,time_index) = repcell(null((H{2,1}*V{1,freq_index(1),time_index(1)})'),[1 LTE_params.IA_freq_granularity LTE_params.IA_time_granularity]);
                    U(3,freq_index,time_index) = repcell(null((H{3,1}*V{1,freq_index(1),time_index(1)})'),[1 LTE_params.IA_freq_granularity LTE_params.IA_time_granularity]);
                end

                channel_error = get_channel_error(H_fft_noisy,ChanMod_output,LTE_params);
            end
    case 'min_WLI'
        P_rec = P;

        switch LTE_params.ChanMod_config.filtering
            case 'BlockFading'
                U_rec = cell(K,LTE_params.Ntot);
                Q = cell(K,LTE_params.Ntot);
                Q_rec = cell(K,LTE_params.Ntot);
        
                for n = 1:LTE_params.Ntot
                    H = get_H(ChanMod_output,n,1,LTE_params);
                    
                    % Initialize precoding matrices
                    for i = 1:K
                        V{i,n} = randn(nTX,d(i));
                    end

                    ii = 1;
                    while(true)
                        for i = 1:K
                            Q{i,n} = zeros(nRX,nRX);
                            for j = 1:K
                                if (j ~= i)
                                    Q{i,n} = Q{i,n} + P(j)/d(j) * H{i,j}*V{j,n}*V{j,n}'*H{i,j}';
                                end
                            end
                        end

                        for i = 1:K
                            [A,B] = eigs(Q{i,n});
                            U{i,n} = A(:,end:-1:end-d(i)+1);

                            p(ii,i) = 0;
                            for k = 1:d(i)
                                p(ii,i) = p(ii,i) + B(end-k+1,end-k+1);
                            end
                            p(ii,i) = p(ii,i) / trace(Q{i,n});
                        end

                        V_rec = U;

                        for j = 1:K
                            Q_rec{j,n} = zeros(nTX,nTX);
                            for i = 1:K
                                if (i ~= j)
                                    Q_rec{j,n} = Q_rec{j,n} + P_rec(i)/d(i) * H{i,j}'*V_rec{i,n}*V_rec{i,n}'*H{i,j}'';
                                end
                            end
                        end

                        for j = 1:K
                            [A,~] = eigs(Q_rec{j,n});
                            U_rec{j,n} = A(:,end:-1:end-d(j)+1);
                        end

                        V = U_rec;

                        WLI(ii) = 0;
                        for k = 1:K
                            for j = 1:K
                                 if (j ~= k)
                                    WLI(ii) = WLI(ii) + P_rec(k)/d(k) * P(j)/d(j) * trace(U{k,n}'*H{k,j}*V{j,n}*V{j,n}'*H{k,j}'*U{k,n});
                                 end
                            end
                        end

                        if (abs(WLI(ii)) < IA_thresh)
                            break;
                        end

                        if (ii >= IA_max_iterations)
                            warning('Maximum number of iterations reached!');
                            break;
                        end

                        ii = ii+1;
                    end
                end
                V = repmat(V,[1 1 LTE_params.Nsub]);
                U = repmat(U,[1 1 LTE_params.Nsub]);
        end

    case 'max_SINR'
        P_rec = P;

        switch LTE_params.ChanMod_config.filtering
            case 'BlockFading'
                U_rec = cell(K,LTE_params.Ntot);
                Q = cell(K,LTE_params.Ntot);
                B = cell(K,max(d),LTE_params.Ntot);
                B_rec = cell(K,max(d),LTE_params.Ntot);
        
                for n = 1:LTE_params.Ntot
                    H = get_H(ChanMod_output,n,1,LTE_params);        
        
                    % Initialize precoding matrices
                    for i = 1:K
                        V{i,n} = randn(nTX,d(i));
                    end

                    ii = 1;
                    while(true)
                        for i = 1:K
                            for k = 1:d(i)
                                B{i,k,n} = zeros(nRX,nRX);
                                for j = 1:K
                                    for dd = 1:d(j)
                                        B{i,k,n} = B{i,k,n} + P(j)/d(j) * H{i,j}*V{j,n}(:,dd)*V{j,n}(:,dd)'*H{i,j}';
                                    end
                                end
                                B{i,k,n} = B{i,k,n} - P(i)/d(i) * H{i,i}*V{i,n}(:,k)*V{i,n}(:,k)'*H{i,i}' + sigma_n2*eye(nRX);
                            end
                        end

                        for i = 1:K
                            Q{i,n} = zeros(nRX,nRX);
                            for j = 1:K
                                if (j ~= i)
                                    Q{i,n} = Q{i,n} + P(j)/d(j) * H{i,j}*V{j,n}*V{j,n}'*H{i,j}';
                                end
                            end
                            [~,D] = eigs(Q{i,n});

                            p(ii,i) = 0;
                            for k = 1:d(i)
                                p(ii,i) = p(ii,i) + D(end-k+1,end-k+1);
                            end
                            p(ii,i) = p(ii,i) / trace(Q{i,n});
                        end    

                        for i = 1:K
                            U{i,n} = zeros(nRX,d(i));
                            for k = 1:d(i)
                                U{i,n}(:,k) = (B{i,k,n}\H{i,i})*V{i,n}(:,k) / norm((B{i,k,n}\H{i,i})*V{i,n}(:,k));
                            end
                        end

                        V_rec = U;

                        for i = 1:K
                            for k = 1:d(i)
                                B_rec{i,k,n} = zeros(nTX,nTX);
                                for j = 1:K
                                    for dd = 1:d(j)
                                        B_rec{i,k,n} = B_rec{i,k,n} + P_rec(j)/d(j) * H{j,i}'*V_rec{j,n}(:,dd)*V_rec{j,n}(:,dd)'*H{j,i}'';
                                    end
                                end
                                B_rec{i,k,n} = B_rec{i,k,n} - P_rec(i)/d(i) * H{i,i}'*V_rec{i,n}(:,k)*V_rec{i,n}(:,k)'*H{i,i}'' + sigma_n2*eye(nTX);
                            end
                        end

                        for i = 1:K
                            U_rec{i,n} = zeros(nTX,d(i));
                            for k = 1:d(i)
                                U_rec{i,n}(:,k) = (B_rec{i,k,n}\H{i,i}')*V_rec{i,n}(:,k) / norm((B_rec{i,k,n}\H{i,i}')*V_rec{i,n}(:,k));
                            end
                        end

                        V = U_rec;

                        SINR(:,:,ii) = zeros(K,d(i));
                        for i = 1:K
                            for k = 1:d(i)
                                for j = 1:K
                                    if (j ~= i)
                                        SINR(i,k,ii) = SINR(i,k,ii) + P(i)/d(i) * ...
                                            U{i,n}(:,k)'*H{i,i}*V{i,n}(:,k)*V{i,n}(:,k)'*H{i,i}'*U{i,n}(:,k) / ...
                                            (U{i,n}(:,k)'*B{i,k,n}*U{i,n}(:,k));
                                    end
                                end
                            end
                        end

                        if (ii > 1) && (norm(10*log10(abs(SINR(:,:,ii)))-10*log10(abs(SINR(:,:,ii-1)))) < IA_thresh)
                            break;
                        end

                        if (ii >= IA_max_iterations)
                            warning('Maximum number of iterations reached!');
                            break;
                        end

                        ii = ii+1;
                    end
                end
                V = repmat(V,[1 1 LTE_params.Nsub]);
                U = repmat(U,[1 1 LTE_params.Nsub]);
        end
        
    case 'SVD_beamforming'
        switch LTE_params.ChanMod_config.filtering
            case 'BlockFading'
                for n = 1:LTE_params.Ntot
                    H = get_H(ChanMod_output,n,1,LTE_params);
                    
                    for i = 1:K
                        [U1,~,V1] = svd(H{i,i});
                        U{i,n} = U1(:,1:d(i));
                        V{i,n} = V1(:,1:d(i));
                    end
                end
                V = repmat(V,[1 1 LTE_params.Nsub]);
                U = repmat(U,[1 1 LTE_params.Nsub]);
        end

    case 'greedy_IA'
        P = ones(K,1);

        switch LTE_params.ChanMod_config.filtering
            case 'BlockFading'
                for n = 1:LTE_params.Ntot
                    H = get_H(ChanMod_output,n,1,LTE_params);
                    
                    for i = 1:K
                        V{i,n} = randn(nTX,d(i));
                    end

                    for ii = 1:10
                        for i = 1:K
                            Q{i,n} = zeros(nRX,nRX);
                            for j = 1:K
                                if (j ~= i)
                                    Q{i,n} = Q{i,n} + P(j)/d(j) * H{i,j}*V{j,n}*V{j,n}'*H{i,j}';
                                end
                            end
                        end

                        for i = 1:K
                            [U1,~,V1] = svd((sigma_n2*eye(nTX) + Q{i,n})^(-1/2) * H{i,i});
                            U{i,n} = U1(:,1:d(i));
                            V{i,n} = V1(:,1:d(i));
                        end
                    end
                end
                V = repmat(V,[1 1 LTE_params.Nsub]);
                U = repmat(U,[1 1 LTE_params.Nsub]);
        end

    case 'random_beamforming'
        switch LTE_params.ChanMod_config.filtering
            case 'BlockFading'
                for n = 1:LTE_params.Ntot
                    for i = 1:K
                        U{i,n} = rand(nRX,d(i)) + 1i*rand(nRX,d(i));
                        U{i,n} = U{i,n} / norm(U{i,n});
                        V{i,n} = rand(nTX,d(i)) + 1i*rand(nTX,d(i));
                        V{i,n} = V{i,n} / norm(V{i,n});
                    end
                end
                V = repmat(V,[1 1 LTE_params.Nsub]);
                U = repmat(U,[1 1 LTE_params.Nsub]);
        end
    
    otherwise
        error('Interference Alignment type not supported.');
end

function H = get_H(H_fft_noisy,n,m,LTE_params)
K = LTE_params.nBS;
nRX = LTE_params.UE_config.nRX;
nTX = LTE_params.BS_config.nTx;
H = cell(K,K);
% for i = 1:K
%     for j = 1:K
%         H{i,j} = reshape(ChanMod_output{j,i}.genie.H_fft(n,m,:,:),nTX,nRX);
%     end
% end
for i = 1:K
    for j = 1:K
        if LTE_params.IA_freq_granularity == 1
            H{i,j} = reshape(H_fft_noisy{j,i}(n,m,:,:),nTX,nRX);
        else
            H{i,j} = reshape(mean(H_fft_noisy{j,i}(n,m,:,:)),nTX,nRX);
%             H{i,j} = reshape(H_fft_noisy{j,i}(n(ceil(length(n)/2)),m,:,:),nTX,nRX);
        end
    end
end

function A_cell = repcell(A,n)
A_cell = cell(n);
if length(n) == 2
    for i = 1:n(1)
        for j = 1:n(2)
            A_cell{i,j} = A;
        end
    end
elseif length(n) == 3
    for i = 1:n(1)
        for j = 1:n(2)
            for k = 1:n(3)
                A_cell{i,j,k} = A;
            end
        end
    end
end

function channel_error = get_channel_error(H_fft_noisy,ChanMod_output,LTE_params)
K = LTE_params.nBS;
channel_error_matrix = zeros(K);
for i = 1:K
    for j = 1:K
        H_fft_IA = zeros(size(ChanMod_output{j,i}.genie.H_fft));
        for n = 1:LTE_params.Ntot/LTE_params.IA_freq_granularity
            freq_index = LTE_params.IA_freq_granularity*(n-1)+1:LTE_params.IA_freq_granularity*n;
            if LTE_params.IA_freq_granularity == 1
                H_fft_IA(freq_index,:,:,:) = repmat(H_fft_noisy{j,i}(freq_index,:,:,:),LTE_params.IA_freq_granularity,1);
            else
                H_fft_IA(freq_index,:,:,:) = repmat(mean(H_fft_noisy{j,i}(freq_index,:,:,:)),LTE_params.IA_freq_granularity,1);
            end
            channel_error_matrix(j,i) = mean(mean(mean(mean(abs(ChanMod_output{j,i}.genie.H_fft).^2)))) ./ mean(mean(mean(mean(abs(H_fft_IA-ChanMod_output{j,i}.genie.H_fft).^2))));
        end
    end
end
channel_error = 10*log10(mean(mean(channel_error_matrix)));