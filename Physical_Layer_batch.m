% File To generate plots from the LTE journal paper
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) 2010 by INTHFT
% www.nt.tuwien.ac.at

clear
clc
tx_mode = 2;


% in order to reprodue the figure from the journal paper set the variables
% N_subframes to 1000
% in case to obtain only course result run it with N_subframes = 10

Simulation_type = 'LTE_journal_paper';
N_subframes = 1000;
N_Ue = 1;
N_Bs = 1;
N_subframes2 = N_subframes;
cqi_i = 1;
channel_type = 'winner_II';
connection_table = true(N_Bs,N_Ue);
Power_diff = 100;

load('channel_2x2_5MHz.mat');
channel = channel{1};
[N_rx N_tx nr_of_delay N_subframe] = size(channel);
SNR_vec = -10:2:30;

LTE_load_parameters;
load('channel_2x2_5MHz.mat');
channel = channel{1};
H_fft = nan(LTE_params.Ntot,N_rx,N_tx,LTE_params.Nsub,N_subframe);


for sub_i = 1:N_subframes
    for rr = 1:N_rx
        for tt = 1:N_tx
            spec = fft([squeeze(channel(rr,tt,:,sub_i)); zeros(LTE_params.Nfft-length(squeeze(channel(rr,tt,:,sub_i))),1)]);
            H_fft(:,rr,tt,:,sub_i) = repmat(spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]),1,LTE_params.Nsub); % remove DC carrier and zeros padded up to size of FFT
        end
    end
end

mi_check = zeros(length(SNR_vec),N_subframes,LTE_params.Nsub);
mi = nan(length(SNR_vec),N_subframes,LTE_params.Nsub);
capacity = nan(length(SNR_vec),N_subframes,LTE_params.Nsub);

for snr_i=1:length(SNR_vec)
    snr_i
    for sub_i=1:N_subframes
        for ss = 1:LTE_params.Nsub
            Hf_temp = H_fft(:,:,:,ss,sub_i);
            sigma_n2 = 10^(-SNR_vec(snr_i)/10);
            [CTxWF, CNoTx, PfTxWF, P, N] = jCapacity(Hf_temp./sqrt(sigma_n2), LTE_params.Ntot);
            
            %             for ii=1:size(Hf_temp,1)
            %                 HH = squeeze(Hf_temp(ii,:,:))./sqrt(sigma_n2);
            %                 mi_check(snr_i,sub_i,ss) = mi_check(snr_i,sub_i,ss) + real(log2(det(eye(size(HH,1))+HH*HH')))/size(Hf_temp,1);
            %             end
            mi(snr_i,sub_i,ss) = CNoTx;
            capacity(snr_i,sub_i,ss) = CTxWF;
        end
    end
end

a = mean(capacity,3);
a = mean(a,2);

b = mean(mi,3);
b = mean(b,2);

%c = mean(mi_check,3);
%c = mean(c,2);

% figure
% plot(SNR_vec,a*5,'r');
% hold on
% plot(SNR_vec,b*5,'b');

%plot(SNR_vec,c*5,'g');
%legend('capacity','mi','mi check')


tx_mode_vec = [2 3];
mutual_mat = nan(length(SNR_vec),length(tx_mode_vec));
for tx_mode_i = 1:length(tx_mode_vec)
    tx_mode = tx_mode_vec(tx_mode_i);
    
    Simulation_type = 'LTE_journal_paper';
    %N_subframes = 1000;
    %N_subframes2 = 1000;
    %cqi_i = 1;
    channel_est = 'PERFECT';
    load('channel_2x2_5MHz.mat');
    channel = channel{1};
    [N_rx N_tx nr_of_delay N_subframe] = size(channel);
    
    LTE_load_parameters;
    load('channel_2x2_5MHz.mat');
    channel = channel{1};
    
    H_fft = nan(LTE_params.Ntot,LTE_params.Nsub,N_rx,N_tx,N_subframe);
    %SNR_vec = -10:10:30;
    %SNR_vec = -10:2:0;
    
    % flat
    % channel = (1/sqrt(2))*(randn(N_rx,N_tx,1,N_subframe) + sqrt(-1)*randn(N_rx,N_tx,1,N_subframe));
    
    %awgn
    % H_temp = [  1,1,1,1;
    %     1,-1,-1,1;
    %     1,-1,1,-1;
    %     1,1,-1,-1];
    % H_temp_small = H_temp(1:N_rx,1:N_tx);
    % channel = repmat(H_temp_small,[1 1 1 N_subframe ]);
    
    
    % tranform of the channel matrix
    
    for sub_i = 1:N_subframe
        for rr = 1:N_rx
            for tt = 1:N_tx
                spec = fft([squeeze(channel(rr,tt,:,sub_i)); zeros(LTE_params.Nfft-length(squeeze(channel(rr,tt,:,sub_i))),1)]);
                H_fft(:,:,rr,tt,sub_i) = repmat(spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]),1,LTE_params.Nsub); % remove DC carrier and zeros padded up to size of FFT
            end
        end
    end
    
    %awgn
    %H_fft = ones(size(H_fft));
    
    
    
    % switch tx_mode
    %     case 1
    %
    %     case 2
    %     case 3
    %     case 4
    %     otherwise
    % end
    
    %find date sub
    %find precoding
    %mutual_vec = zeros(N_subframe,length(SNR_vec));
    mutual_vec = zeros(N_subframes2,length(SNR_vec));
    
    for snr_i = 1:length(SNR_vec)
        snr_i
        %sigma_n2 = 10^(-SNR_vec(snr_i)/10);
        sigma_n2 = N_rx*10^(-SNR_vec(snr_i)/10);
        for sub_i = 1:N_subframes2
            
            y_tx_assembled = nan(LTE_params.Ntot,LTE_params.Nsub);
            
            sub_mod = mod(sub_i,10);
            if sub_mod==0
                sub_mod = 10;
            end
            
            RefSym          = LTE_params.Reference_Signal(1,sub_mod).RefSym;
            RefMapping      = LTE_params.Reference_Signal(1,sub_mod).RefMapping;
            total_no_refsym = LTE_params.Reference_Signal(1,sub_mod).total_no_refsym;     % total number of reference symbols per slot per resource block
            NoData_indices  = LTE_params.Reference_Signal(1,sub_mod).NoData_indices;
            
            PrimSync         = LTE_params.Sync_Signal(1,sub_mod).PrimSync;
            PrimMapping      = LTE_params.Sync_Signal(1,sub_mod).PrimMapping;
            SecSync          = LTE_params.Sync_Signal(1,sub_mod).SecSync;
            SecMapping       = LTE_params.Sync_Signal(1,sub_mod).SecMapping;
            SyncUsedElements = LTE_params.Sync_Signal(1,sub_mod).SyncUsedElements;
            ResMapping = LTE_params.Sync_Signal(1,sub_mod).ResMapping;    % reserved REs, nothing is transmitted here!
            
            %CHusedElements = CHusedElements + LTE_params.Sync_Signal(1,subframe_corr).ResUsedElements;
            %     CHnsyms_subframe = CHnsyms_subframe + LTE_params.Sync_Signal(1,subframe_corr).ResNElements;
            
            
            CHmapping = zeros(size(NoData_indices));
            CHusedElements = zeros(LTE_params.Nrb*2,1);
            CHusedElements = CHusedElements + LTE_params.Sync_Signal(1,sub_mod).ResUsedElements;
            CHnsyms_subframe = 0;
            CHmapping = CHmapping | ResMapping;
            %     LTE_params.usePBCH = 1;
            %     LTE_params.usePDCCH = 1;
            
            if LTE_params.usePBCH
                CHmapping = CHmapping | LTE_params.PBCH(1,sub_mod).Mapping;
                CHusedElements = CHusedElements + LTE_params.PBCH(1,sub_mod).UsedElements;
                %     CHnsyms_subframe = CHnsyms_subframe + LTE_params.PBCH(1,subframe_corr).N_Elements;
            end
            if LTE_params.usePDCCH
                CHmapping = CHmapping | LTE_params.PDCCH(1,sub_mod).Mapping;
                CHusedElements = CHusedElements + LTE_params.PDCCH(1,sub_mod).UsedElements;
                %     CHnsyms_subframe = CHnsyms_subframe + LTE_params.PDCCH(1,subframe_corr).N_Elements;
            end
            
            
            
            no_data = LTE_params.Reference_Signal(sub_mod).NoData_indices | LTE_params.Sync_Signal(sub_mod).PrimMapping | LTE_params.Sync_Signal(sub_mod).SecMapping;
            total_no_refsym = LTE_params.Reference_Signal(1,sub_mod).total_no_refsym;
            
            rb_numbers = 1:25*2; % 300subcarriers/12 * 2
            if(sub_mod == 1 || sub_mod == 6)
                rb_tx_symbols = LTE_params.Nsc*LTE_params.Ns - total_no_refsym - SyncUsedElements(true(25,2)) - CHusedElements(true(25,2));
                index = 0;
                for ii = 1:length(rb_numbers)
                    if(rb_numbers(ii) > LTE_params.Nrb)
                        y_tx_assembled_temp = zeros(LTE_params.Nsc,LTE_params.Ns);
                        y_tx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1-LTE_params.Nrb)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,8:8+LTE_params.Ns-1))) = [index+1:index+rb_tx_symbols(ii)];
                        y_tx_assembled((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end) = y_tx_assembled_temp;
                    else
                        y_tx_assembled_temp = zeros(LTE_params.Nsc,LTE_params.Ns);
                        y_tx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)...
                            + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns))) = [index+1:index+rb_tx_symbols(ii)];
                        y_tx_assembled((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) = y_tx_assembled_temp;
                    end
                    index = index + rb_tx_symbols(ii);
                end
            else
                rb_tx_symbols = LTE_params.Nsc*LTE_params.Ns - total_no_refsym;% - CHusedElements(UE_mapping);
                index = 0;
                for ii = 1:length(rb_numbers)
                    if(rb_numbers(ii) > LTE_params.Nrb)
                        y_tx_assembled_temp = zeros(LTE_params.Nsc,LTE_params.Ns);
                        y_tx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))) = [index+1:index+rb_tx_symbols];
                        y_tx_assembled((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end) = y_tx_assembled_temp;
                    else
                        y_tx_assembled_temp = zeros(LTE_params.Nsc,LTE_params.Ns);
                        y_tx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))) = [index+1:index+rb_tx_symbols];
                        y_tx_assembled((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) = y_tx_assembled_temp;
                    end
                    index = index + rb_tx_symbols;
                end
            end
            
            
            switch tx_mode
                case 1
                    mapping = logical(y_tx_assembled);
                    H_fft_temp = squeeze(H_fft(:,:,1,1,sub_i));
                    channel_used = H_fft_temp(mapping);
                    mutal_info_temp = 0;
                    for channel_i = 1:length(channel_used)
                        mutal_info_temp = mutal_info_temp + log2(det(eye(N_rx)+(1/sigma_n2)*channel_used(channel_i)*conj(channel_used(channel_i))));
                    end
                    mutal_info_temp = mutal_info_temp*LTE_params.SubcarrierSpacing;
                    mutual_vec(sub_i,snr_i) = mutal_info_temp;
                case 2
                    mapping = logical(y_tx_assembled);
                    channel_used = nan(sum(sum(mapping)),N_rx,N_tx);
                    for rr = 1:N_rx
                        for tt = 1:N_tx
                            H_fft_temp = squeeze(H_fft(:,:,rr,tt,sub_i));
                            channel_used(:,rr,tt) = H_fft_temp(mapping);
                        end
                    end
                    
                    if N_tx==2
                        mutal_info_temp = 0;
                        for channel_i = 1:length(channel_used)
                            h_1 = channel_used(channel_i,:,1).';
                            h_2 = channel_used(channel_i,:,2).';
                            virtual_channel_temp = [h_1 -h_2; conj(h_2) conj(h_1) ];
                            mutal_info_temp = mutal_info_temp + log2(det(eye(N_rx*2)+(1/sigma_n2)*virtual_channel_temp*(virtual_channel_temp')));
                        end
                        mutal_info_temp = mutal_info_temp*LTE_params.SubcarrierSpacing/2;
                        mutual_vec(sub_i,snr_i) = mutal_info_temp;
                    elseif N_tx==4
                        mutal_info_temp = 0;
                        for channel_i = 1:length(channel_used)
                            h_1 = channel_used(channel_i,:,1).';
                            h_2 = channel_used(channel_i,:,2).';
                            h_3 = channel_used(channel_i,:,3).';
                            h_4 = channel_used(channel_i,:,4).';
                            virtual_channel_temp = [h_1 -h_2 zeros(N_rx,1) zeros(N_rx,1); conj(h_2) conj(h_1) zeros(N_rx,1) zeros(N_rx,1);zeros(N_rx,1) zeros(N_rx,1) h_3 -h_4; zeros(N_rx,1) zeros(N_rx,1) conj(h_4) conj(h_3)];
                            mutal_info_temp = mutal_info_temp + log2(det(eye(N_rx*4)+(1/sigma_n2)*virtual_channel_temp*(virtual_channel_temp')));
                        end
                        mutal_info_temp = mutal_info_temp*LTE_params.SubcarrierSpacing/4;
                        mutual_vec(sub_i,snr_i) = mutal_info_temp;
                    end
                case 3
                    mapping = logical(y_tx_assembled);
                    channel_used = nan(sum(sum(mapping)),N_rx,N_tx);
                    for rr = 1:N_rx
                        for tt = 1:N_tx
                            H_fft_temp = squeeze(H_fft(:,:,rr,tt,sub_i));
                            channel_used(:,rr,tt) = H_fft_temp(mapping);
                        end
                    end
                    %precoding
                    if N_tx==2
                        [Z W U D] = LTE_common_get_precoding_matrix(tx_mode,N_tx,1,min(N_rx,N_tx),LTE_params);
                        PRE = zeros(N_tx,min(N_rx,N_tx),2);
                        PRE(:,:,1) = W*D*U;
                        PRE(:,:,2) = W*D^2*U;
                        
                        effective_channel = nan(sum(sum(mapping)),min(N_rx,N_tx),min(N_rx,N_tx));
                        mutal_info_temp = 0;
                        
                        channel_used_shifted = shiftdim(channel_used,1);
                        for channel_i = 1:size(channel_used,1)
                            if mod(channel_i,2)==1
                                pre_used = PRE(:,:,1);
                            else
                                pre_used = PRE(:,:,2);
                            end
                            
                            channel_temp = channel_used_shifted(:,:,channel_i)*pre_used;
                            effective_channel(channel_i,:,:) = channel_temp;
                            mutal_info_temp = mutal_info_temp + log2(det(eye(N_rx)+(1/sigma_n2)*channel_temp*(channel_temp')));
                            
                        end
                        mutal_info_temp = mutal_info_temp*LTE_params.SubcarrierSpacing;
                        mutual_vec(sub_i,snr_i) = mutal_info_temp;
                        
                    elseif N_tx==4
                        [Z W U D] = LTE_common_get_precoding_matrix(tx_mode,N_tx,[12,13,14,15],min(N_rx,N_tx),LTE_params);
                        
                        l = 1:length(channel_used);
                        k = mod(floor((l-1)/min(N_rx,N_tx)),4)+1;
                        p = mod(l-1,min(N_rx,N_tx))+1;
                        period = lcm(seqperiod(k),seqperiod(p));
                        PRE = zeros(4,min(N_rx,N_tx),period);
                        for i = 1:period
                            PRE(:,:,i)= W(:,:,k(i))*D^(p(i))*U;
                        end
                        
                        effective_channel = nan(sum(sum(mapping)),min(N_rx,N_tx),min(N_rx,N_tx));
                        mutal_info_temp = 0;
                        
                        channel_used_shifted = shiftdim(channel_used,1);
                        for channel_i = 1:size(channel_used,1)
                            used_index = mod(channel_i,period);
                            if used_index==0
                                used_index = period;
                            end
                            pre_used = PRE(:,:,used_index);
                            channel_temp = channel_used_shifted(:,:,channel_i)*pre_used;
                            effective_channel(channel_i,:,:) = channel_temp;
                            mutal_info_temp = mutal_info_temp + log2(det(eye(N_rx)+(1/sigma_n2)*channel_temp*channel_temp'));
                            
                        end
                        mutal_info_temp = mutal_info_temp*LTE_params.SubcarrierSpacing;
                        mutual_vec(sub_i,snr_i) = mutal_info_temp;
                        
                    end
                    
                case 4
                otherwise
            end
            
            
        end
    end
    mean(mutual_vec,1);
    mutual_vec = mutual_vec*14/15; %cyclic prefix factor
    mutual_vec = mutual_vec/14;
    
    mutual_mat(:,tx_mode_i) = mean(mutual_vec);
    %figure
    %plot(SNR_vec,mean(mutual_vec)/1e6,'r')
end




%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;

%% SNR setting
% SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
% SNR_stepsize = 1;
% SNR_window = 0.25;
% power = [];
% noise = [];
Simulation_type = 'LTE_journal_paper';     %'SUSISO'
                                    %'MUSISO'
                                    %'SUMIMO'
                                    %'MUMIMO'
                                    %'SUSISO_quick_test'
                                    %'SUSISO_BLER_curves_batch'
                                    %'SUSISO_best_cqi'
                                    %'SUMIMO_quick_test'
                                    %'winner_model_example'
                                    %'wsa_2010_michal'
                                    
                                    
%  counti = 1;
%  channel_estimation_error_freq_depend = zeros(72,14,1,2,500);
% Hsave = zeros(72,1000,32);
%N_subframes = 10;

%SNR_vec = -10:1:30;
%tx_mode_vec = [2 3];
%N_rx = 2; N_tx = 2;
tho_matrix = nan(length(tx_mode_vec),N_subframes,length(SNR_vec),15);
channel_est = 'PERFECT';
equalizer = 'SSD';
%% Actual simulations
for tx_mode_i = 1:length(tx_mode_vec)
    tx_mode = tx_mode_vec(tx_mode_i);
    for cqi_i = 1:15


          LTE_load_parameters;  % Single User Multiple Input Multiple Output

        LTE_sim_main
        tho_matrix(tx_mode_i,:,:,cqi_i) = sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6;

    end
end

%% plot
figure
 
plot(SNR_vec,a*5,'color','m','Linestyle','-','Marker','.'); %multiply by 5, due to the 5MHz bandwidth
hold on
plot(SNR_vec,b*5,'color','black','Linestyle','-','Marker','.');

plot(SNR_vec,mutual_mat(:,2)/1e6,'color','b','Linestyle','-.','Marker','.');
plot(SNR_vec,mutual_mat(:,1)/1e6,'color','r','Linestyle','-.','Marker','.');

tho = squeeze(max(tho_matrix,[],4)); %selection of the optimal cqi
tho_mean = squeeze(mean(tho,2));

plot(SNR_vec,tho_mean(2,:),'color','b','Linestyle',':','Marker','.');
plot(SNR_vec,tho_mean(1,:),'color','r','Linestyle',':','Marker','.');

xlabel('SNR [dB]')
ylabel('Throughput [Mbit/s]')
grid on
legend('Capacity','MI','AMI OLSM','AMI TxD','OSLM','TxD','Location','NorthWest')

