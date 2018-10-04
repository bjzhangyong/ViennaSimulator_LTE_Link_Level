% Reproducibility batch script of paper "Interference Alignment in UMTS
% Long Term Evolution"
% Author: Joerg Reitterer, jreitter@nt.tuwien.ac.at
% (c) 2011 by INTHFT
% www.nt.tuwien.ac.at

clear all;
close all;

LineWidth = 1;
MarkerSize = 4;

%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;


%% Figure 5(a)
Simulation_type = 'IA'; 
N_Ue = 1;
N_Bs = 3;
tx_mode = 4;
N_rx = 2;
N_tx = 2;
channel_type = 'VehA';
SNR_vec = -10:2.5:30;
SIR_vec = 0:2.5:30;
connection_table = logical(eye(3));
IA_type = 'closed_form'; % Type of IA; closed_form, min_WLI, or max_SINR
IA_streams = [2 2 2]; % Number of spatial streams
IA_thresh = 1e-5; % IA threshold for stopping iterative algorithms
IA_max_iterations = 1000;
IA_freq_granularity_vec = 1; % Spectral IA granularity;
                             % 1 ... alignment is calclulated for each subcarrier
                             % 12 ... alignment is calclulated for each resource block
                             % etc.
IA_time_granularity = 14; % IA time granularity;
                          % 1 ... alignment is calclulated for each OFDM symbol
                          % 7 ... alignment is calculated for each slot
                          % 14 ... alignment is calcualted for each subframe
                          % etc.
receiver = 'ZF';
channel_estimation_method = 'PERFECT';
IA_sigma_H2_E2_ratio_vec = Inf; % Channel measurement error;
                                % Inf ... perfect channel knowledge
scheduler_type = 'round robin';
scheduler_assignment = 'static';
user_speed_vec = 0;
filtering = 'BlockFading';

cell_throughput_useful = zeros(length(SNR_vec),length(SIR_vec));
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
                cell_throughput_useful(:,SIR_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
            end
        end       
    end
end

figure;
surf(SNR_vec(1,:),SIR_vec,cell_throughput_useful.');
view(-115,30)
xlim([-10 30]);
zlim([0 5]);
xlabel('SNR [dB]');
ylabel('SIR [dB]');
zlabel('Throughput [Mbit/s]');
title_text = ['LTE CLSM, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', d_i = [2, 2, 2], CQI = ' num2str(cqi_i) ', ' channel_type];
title(title_text);
grid on;
box on;
caxis_ref = caxis;


%% Figure 5(b)
Simulation_type = 'IA';
tx_mode = 6;
IA_streams = [1 1 1];
SNR_vec = -10:2.5:30;
SIR_vec = 0:2.5:30;

cell_throughput_useful = zeros(length(SNR_vec),length(SIR_vec));
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
                cell_throughput_useful(:,SIR_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
            end
        end       
    end
end

figure;
surf(SNR_vec(1,:),SIR_vec,cell_throughput_useful.');
view(-115,30)
xlim([-10 30]);
zlim([0 5]);
xlabel('SNR [dB]');
ylabel('SIR [dB]');
zlabel('Throughput [Mbit/s]');
title_text = ['LTE IA, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', d_i = [1, 1, 1], CQI = ' num2str(cqi_i) ', ' channel_type];
title(title_text);
grid on;
box on;
caxis(caxis_ref);


% %% Figure 6(a)
SNR_vec = 30;
SIR_vec = 0;
IA_freq_granularity_vec = [1 2 3 4 6 12 24 36 72];
channel_type = 'flat Rayleigh';

cell_throughput_useful_flat_Rayleigh = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 15
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_flat_Rayleigh(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

channel_type = 'PedA';

cell_throughput_useful_PedA = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 15
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_PedA(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

channel_type = 'VehA';

cell_throughput_useful_VehA = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 15
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_VehA(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

figure;
hold on;
plot(IA_freq_granularity_vec,cell_throughput_useful_flat_Rayleigh/cell_throughput_useful_flat_Rayleigh(1)*100,'bo-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
plot(IA_freq_granularity_vec,cell_throughput_useful_PedA/cell_throughput_useful_PedA(1)*100,'rs-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
plot(IA_freq_granularity_vec,cell_throughput_useful_VehA/cell_throughput_useful_VehA(1)*100,'gd-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
hold off;
grid on;
box on;
ylim([0 105]);
title_text = ['LTE IA, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', SNR = ' num2str(SNR_vec(1)) ' dB, SIR = ' num2str(SIR_vec(1)) ' dB, CQI = ' num2str(cqi_i(1))];
title(title_text);
xlabel('\xi_f')
ylabel('Throughput [%]')
legend('Flat Rayleigh','ITU-T PedA','ITU-T VehA','Location','SouthEast');


%% Figure 6(b)
SNR_vec = 15;
SIR_vec = 0;
IA_freq_granularity_vec = [1 2 3 4 6 12 24 36 72];
channel_type = 'flat Rayleigh';

cell_throughput_useful_flat_Rayleigh = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_flat_Rayleigh(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

channel_type = 'PedA';

cell_throughput_useful_PedA = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_PedA(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

channel_type = 'VehA';

cell_throughput_useful_VehA = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_VehA(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

figure;
hold on;
plot(IA_freq_granularity_vec,cell_throughput_useful_flat_Rayleigh/cell_throughput_useful_flat_Rayleigh(1)*100,'bo-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
plot(IA_freq_granularity_vec,cell_throughput_useful_PedA/cell_throughput_useful_PedA(1)*100,'rs-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
plot(IA_freq_granularity_vec,cell_throughput_useful_VehA/cell_throughput_useful_VehA(1)*100,'gd-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
hold off;
grid on;
box on;
ylim([0 105]);
title_text = ['LTE IA, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', SNR = ' num2str(SNR_vec(1)) ' dB, SIR = ' num2str(SIR_vec(1)) ' dB, CQI = ' num2str(cqi_i(1))];
title(title_text);
xlabel('\xi_f')
ylabel('Throughput [%]')
legend('Flat Rayleigh','ITU-T PedA','ITU-T VehA','Location','East');


%% Figure 6(c)
SNR_vec = 0;
SIR_vec = 0;
IA_freq_granularity_vec = [1 2 3 4 6 12 24 36 72];
channel_type = 'flat Rayleigh';

cell_throughput_useful_flat_Rayleigh = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 3
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_flat_Rayleigh(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

channel_type = 'PedA';

cell_throughput_useful_PedA = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 3
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_PedA(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

channel_type = 'VehA';

cell_throughput_useful_VehA = zeros(length(IA_freq_granularity_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 3
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful_VehA(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

figure;
hold on;
plot(IA_freq_granularity_vec,cell_throughput_useful_flat_Rayleigh/cell_throughput_useful_flat_Rayleigh(1)*100,'bo-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
plot(IA_freq_granularity_vec,cell_throughput_useful_PedA/cell_throughput_useful_PedA(1)*100,'rs-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
plot(IA_freq_granularity_vec,cell_throughput_useful_VehA/cell_throughput_useful_VehA(1)*100,'gd-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
hold off;
grid on;
box on;
ylim([0 105]);
title_text = ['LTE IA, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', SNR = ' num2str(SNR_vec(1)) ' dB, SIR = ' num2str(SIR_vec(1)) ' dB, CQI = ' num2str(cqi_i(1))];
title(title_text);
xlabel('\xi_f')
ylabel('Throughput [%]')
legend('Flat Rayleigh','ITU-T PedA','ITU-T VehA','Location','SouthEast');


%% Figure 7
SNR_vec = 15;
SIR_vec = [0:2.5:10 20 30];
IA_freq_granularity_vec = [1 2 3 4 6 12 24 36 72];
channel_type = 'VehA';

cell_throughput_useful = zeros(length(SIR_vec),length(IA_freq_granularity_vec));
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
                cell_throughput_useful(SIR_i,IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
            end
        end
        
    end
end

figure;
plot(IA_freq_granularity_vec,cell_throughput_useful/cell_throughput_useful(1,1)*100,'o-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
grid on;
box on;
ylim([0 105]);
title_text = ['LTE IA, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', SNR = ' num2str(SNR_vec(1)) ' dB, CQI = ' num2str(cqi_i(1)) ', ' channel_type];
title(title_text);
xlabel('\xi_f')
ylabel('Throughput [%]')
for k = 1:length(SIR_vec)
    legend_labels{k} = ['SIR = ' num2str(SIR_vec(k)) ' dB'];
end
hleg = legend(legend_labels,'Location','SouthWest');


%% Figure 8
SNR_vec = 15;
SIR_vec = 0;
IA_freq_granularity_vec = 1;
IA_time_granularity = 14;

IA_sigma_H2_E2_ratio_vec = [0:2.5:30 Inf];

channel_type = 'flat Rayleigh';

cell_throughput_useful = zeros(length(IA_sigma_H2_E2_ratio_vec),1);
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
            cell_throughput_useful(IA_sigma_H2_E2_ratio_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
        end
        
    end
end

hfig = figure;
hold on;
plot(IA_sigma_H2_E2_ratio_vec,ones(size(IA_sigma_H2_E2_ratio_vec))*cell_throughput_useful(end)/cell_throughput_useful(end)*100,'k--','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
plot(IA_sigma_H2_E2_ratio_vec,cell_throughput_useful/cell_throughput_useful(end)*100,'bo-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
hold off;
grid on;
box on;
ylim([0 105]);
title_text = ['LTE IA, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', SNR = ' num2str(SNR_vec(1)) ' dB, CQI = ' num2str(cqi_i(1)) ', Flat Rayleigh'];
title(title_text);
xlabel('\sigma_H^2/\sigma_E^2 [dB]')
ylabel('Throughput [%]')
hleg = legend('Perfect CSI','Noisy CSI','Location','SouthEast');


%% Figure 9
SNR_vec = 15;
SIR_vec = 0;
IA_freq_granularity_vec = 1;
IA_freq_granularity = IA_freq_granularity_vec;
IA_time_granularity_vec = [1 7 14];

IA_sigma_H2_E2_ratio_vec = Inf;

channel_type = 'VehA';
filtering = 'FastFading';
user_speed_vec = (0:50:50)./3.6;

cell_throughput_useful = zeros(length(IA_time_granularity_vec),length(user_speed_vec));
% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_time_granularity_i = 1:length(IA_time_granularity_vec)
        IA_time_granularity = IA_time_granularity_vec(IA_time_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 500;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                end
            end
        end
        cell_throughput_useful(IA_time_granularity_i,user_speed_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
end

hfig = figure;
plot(user_speed_vec*3.6,cell_throughput_useful/cell_throughput_useful(1,1)*100,'o-','LineWidth',LineWidth,'MarkerSize',MarkerSize,'MarkerFaceColor','w');
grid on;
box on;
xlim([0 75]);
ylim([0 105]);
set(gca,'XTick',0:25:75);
title_text = ['LTE IA, K = ' num2str(N_Bs) ', ' num2str(N_tx) 'x' num2str(N_rx) ', SNR = ' num2str(SNR_vec(1)) ' dB, SIR = ' num2str(SIR_vec(1)) ' dB, CQI = ' num2str(cqi_i(1)) ', ' channel_type];
title(title_text);
xlabel('User Velocity v [km/h]')
ylabel('Throughput [%]')
hleg = legend('\xi_t = 1','\xi_t = 7','\xi_t =  14','Location','SouthWest');