% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear
%clear global
%close all
%clc

%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;

%% SNR setting
% SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
% SNR_stepsize = 1;
% SNR_window = 0.25;
% power = [];
% noise = [];
Simulation_type = 'IA';     %'SUSISO'
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

N_Ue = 1;
N_Bs = 3;
tx_mode = 6;
N_rx = 2;
N_tx = 2;
channel_type = 'VehA';
% SNR_vec = -10:2.5:30;
SNR_vec = 15;
SIR_vec = 0;
connection_table = logical(eye(3));
IA_type = 'closed_form'; % Type of IA; closed_form, min_WLI, or max_SINR
IA_streams = [1 1 1]; % Number of spatial streams
IA_thresh = 1e-5; % IA threshold for stopping iterative algorithms
IA_max_iterations = 1000;
%IA_freq_granularity_vec = [1 2 3 4 6 12 24 36 72];
%IA_freq_granularity_vec = [1 12 36 72];
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

user_speed_vec = (0:10:50)./3.6;
% user_speed_vec = 30/3.6;
filtering = 'FastFading';

%% Actual simulations
for user_speed_i = 1:length(user_speed_vec)
    user_speed = user_speed_vec(user_speed_i);
    for IA_freq_granularity_i = 1:length(IA_freq_granularity_vec)
        IA_freq_granularity = IA_freq_granularity_vec(IA_freq_granularity_i);
        for IA_sigma_H2_E2_ratio_i = 1:length(IA_sigma_H2_E2_ratio_vec)
            IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio_vec(IA_sigma_H2_E2_ratio_i);
            for SIR_i = 1:length(SIR_vec)
                Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relative signal power from one eNodeB to the other 
                for cqi_i = 9
                    N_subframes = 10;
                    LTE_load_parameters;  % Single User Multiple Input Multiple Output
                    LTE_params.max_HARQ_retransmissions = 0;
                    LTE_params.channel_noise = 0;
                    LTE_params.ICI_est_type = 'PERFECT';
                    LTE_params.show_plots = false;
                    LTE_sim_main
                    % Code to generate the output filename
        %             output_filename = ['LTE_IA_' num2str(tx_mode) '_' num2str(N_Bs) '_' num2str(N_tx) 'x' num2str(N_rx) '_' num2str(IA_streams, '%i') '_' IA_type '_' num2str(SIR_vec(SIR_i)) '_' num2str(cqi_i) '_' num2str(N_subframes) '_' channel_type '_' LTE_params.UE_config.receiver]
            %         save(fullfile([output_filename '.mat']));
                end
            end
    %         throughput(IA_sigma_H2_E2_ratio_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
        end
%         IA_channel_error_vec(IA_freq_granularity_i) = IA_channel_error;
%         throughput(IA_freq_granularity_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
    end
    throughput(user_speed_i) = mean(sum(simulation_results.cell_specific.throughput_useful,3)/LTE_params.Tsubframe/1e6);
end

% plot(IA_sigma_H2_E2_ratio_vec,throughput,'bo-')
% plot(IA_freq_granularity_vec,throughput,'bo-');
% grid on;

throughput
plot(user_speed_vec*3.6,throughput,'bo-');
grid on;

% N_Ue = 1;
% N_Bs = 3;
% tx_mode = 2;
% N_rx = 2;
% N_tx = 2;
% channel_type = 'PedA';
% SNR_vec = -10:2.5:30;
% SIR_vec = Inf;
% connection_table = logical(eye(3));
% IA_type = 'closed_form';
% IA_streams = [1 1 1];
% 
% %% Actual simulations
% for SIR_i = 1:length(SIR_vec)
%     Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relativ signal power from one eNodeB to the other 
%     for cqi_i = 1:15
%         N_subframes = 1000;
%         LTE_load_parameters;  % Single User Multiple Input Multiple Output
%         LTE_params.show_plots = false;
%         LTE_sim_main
%         % Code to generate the output filename
%         output_filename = ['LTE_IA_' num2str(tx_mode) '_' num2str(N_Bs) '_' num2str(N_tx) 'x' num2str(N_rx) '_' num2str(IA_streams, '%i') '_' IA_type '_' num2str(SIR_vec(SIR_i)) '_' num2str(cqi_i) '_' num2str(N_subframes) '_' channel_type]
%         save(fullfile([output_filename '.mat']));
%     end
% end
% 
% N_Ue = 1;
% N_Bs = 3;
% tx_mode = 3;
% N_rx = 2;
% N_tx = 2;
% channel_type = 'PedA';
% SNR_vec = -10:2.5:30;
% SIR_vec = Inf;
% connection_table = logical(eye(3));
% IA_type = 'closed_form';
% IA_streams = [1 1 1];
% 
% %% Actual simulations
% for SIR_i = 1:length(SIR_vec)
%     Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relativ signal power from one eNodeB to the other 
%     for cqi_i = 1:15
%         N_subframes = 1000;
%         LTE_load_parameters;  % Single User Multiple Input Multiple Output
%         LTE_params.show_plots = false;
%         LTE_sim_main
%         % Code to generate the output filename
%         output_filename = ['LTE_IA_' num2str(tx_mode) '_' num2str(N_Bs) '_' num2str(N_tx) 'x' num2str(N_rx) '_' num2str(IA_streams, '%i') '_' IA_type '_' num2str(SIR_vec(SIR_i)) '_' num2str(cqi_i) '_' num2str(N_subframes) '_' channel_type]
%         save(fullfile([output_filename '.mat']));
%     end
% end
% 
% N_Ue = 1;
% N_Bs = 3;
% tx_mode = 4;
% N_rx = 2;
% N_tx = 2;
% channel_type = 'PedA';
% SNR_vec = -10:2.5:30;
% SIR_vec = Inf;
% connection_table = logical(eye(3));
% IA_type = 'closed_form';
% IA_streams = [1 1 1];
% 
% %% Actual simulations
% for SIR_i = 1:length(SIR_vec)
%     Power_diff = SIR_vec(SIR_i)*ones(N_Bs,1); % difference of the relativ signal power from one eNodeB to the other 
%     for cqi_i = 1:15
%         N_subframes = 1000;
%         LTE_load_parameters;  % Single User Multiple Input Multiple Output
%         LTE_params.show_plots = false;
%         LTE_sim_main
%         % Code to generate the output filename
%         output_filename = ['LTE_IA_' num2str(tx_mode) '_' num2str(N_Bs) '_' num2str(N_tx) 'x' num2str(N_rx) '_' num2str(IA_streams, '%i') '_' IA_type '_' num2str(SIR_vec(SIR_i)) '_' num2str(cqi_i) '_' num2str(N_subframes) '_' channel_type]
%         save(fullfile([output_filename '.mat']));
%     end
% end

% shutdown(10)
% simulation_results.cell_specific.BER_uncoded_overall
% simulation_results.cell_specific.BER_coded_overall
% mean(sum(simulation_results.cell_specific.throughput_uncoded,3)/LTE_params.Tsubframe/1e6)
% mean(sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6)

% confidence_intervals_color =  [0.8 0.8 0.8];

% %% Cell throughtput plot
% cell_throughput_plot_figure = figure;
% 
% % Plot total throughput (sum of all streams)
% cell_throughput_coded = sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6;
% plot(SNR_vec,mean(cell_throughput_coded,1),'.-b','Markersize',5);
% hold on
% cell_throughput_uncoded = sum(simulation_results.cell_specific.throughput_uncoded,3)/LTE_params.Tsubframe/1e6;
% plot(SNR_vec,mean(cell_throughput_uncoded,1),'.-r','Markersize',5);
% % 
% % LTE_plot_confidence_intervals(SNR_vec, confidence_intervals_color ,cell_throughput_plot_figure,'mean', cell_throughput_coded);
% % LTE_plot_confidence_intervals(SNR_vec, confidence_intervals_color ,cell_throughput_plot_figure,'mean', cell_throughput_uncoded);
% 
% legend('cell coded throughput','cell uncoded throughput','Location','best');
% xlabel('SNR [dB]');
% ylabel('Throughput [Mbit/s]');
% title('Cell throughput');
% hold off
% grid on