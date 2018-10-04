% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear
clear global
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
Simulation_type = 'vtc_2011_spring_michal';     %'SUSISO'
%'MUSISO'
%'SUMIMO'
%'MUMIMO'
%'SUSISO_quick_test'
%'SUSISO_BLER_curves_batch'
%'SUSISO_best_cqi'
%'SUMIMO_quick_test'
%'winner_model_example'
%'wsa_2010_michal'



SNR_vec = 30;
speed_vec = [0:25:350]/3.6;
channel_noise_vec = [0];
subcarrier_intervals_vec = [5];
order_vec = [1:5];
estimator_vec = 1:5;
N_subframes = 1000;
tho = nan(4,length(estimator_vec),length(speed_vec),length(channel_noise_vec),length(subcarrier_intervals_vec),length(order_vec));
mse = nan(4,length(estimator_vec),length(speed_vec),length(channel_noise_vec),length(subcarrier_intervals_vec),length(order_vec));
post_sinr = nan(4,length(estimator_vec),length(speed_vec),length(channel_noise_vec),length(subcarrier_intervals_vec),length(order_vec));


signal = nan(length(estimator_vec),length(speed_vec),length(channel_noise_vec),length(subcarrier_intervals_vec),length(order_vec),N_subframes);
interference = nan(length(estimator_vec),length(speed_vec),length(channel_noise_vec),length(subcarrier_intervals_vec),length(order_vec),N_subframes);
SINR_optimal = nan(length(estimator_vec),length(speed_vec),length(channel_noise_vec),length(subcarrier_intervals_vec),length(order_vec),N_subframes);
%mean_channel_matrix = nan(length(speed_vec),N_subframes,72,72);

%LTE_params.ICI_est_type = '1st';%1st_ort, thomas_1st,1st

%% Actual simulations
for ici_est_i = 1:4
    switch ici_est_i
        case 1
            LTE_params.ICI_est_type = 'PERFECT';%1st_ort, thomas_1st,1st
        case 2
            LTE_params.ICI_est_type = '1st';%1st_ort, thomas_1st,1st
        case 3
            LTE_params.ICI_est_type = '1st_ort';%1st_ort, thomas_1st,1st
        case 4
            LTE_params.ICI_est_type = 'thomas_1st';%1st_ort, thomas_1st,1st
    end
         
    for est_i = 2:2
        switch est_i
            case 1
                estimator = 'LS_michal'
            case 2
                estimator = 'MMSE_Rayleigh2'
            case 3
                estimator = 'PERFECT'
            case 4
                estimator = 'MMSE'
            case 5
                estimator = 'LS'
        end
        for sub_i = 1:length(subcarrier_intervals_vec)
            subcarrier_intervals = subcarrier_intervals_vec(sub_i)
            for channel_noise_i = 1:length(channel_noise_vec)
                for order_i = 1:length(order_vec)
                    for speed_i = 1:length(speed_vec)
                        for cqi_i = 14
                            
                            speed = speed_vec(speed_i)
                            %N_subframes = 5;
                            
                            
                            LTE_load_parameters  % Multi User Multiple Input Multiple Output
                            LTE_params.channel_noise = channel_noise_vec(channel_noise_i);
                            LTE_params.subcarrier_intervals = subcarrier_intervals;
                            LTE_params.order = order_vec(order_i);
                            
                            LTE_sim_main;
                            
                            tho(ici_est_i,est_i,speed_i,channel_noise_i,sub_i,order_i) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
                            mse(ici_est_i,est_i,speed_i,channel_noise_i,sub_i,order_i) = simulation_results.cell_specific.MSE_overall;
                            post_sinr(ici_est_i,est_i,speed_i,channel_noise_i,sub_i,order_i) = mean(simulation_results.cell_specific.PE_SINR_overall);
                            %                         signal(est_i,speed_i,channel_noise_i,sub_i,order_i,:) = signal_power;
                            %                         interference(est_i,sp
                            %                         eed_i,channel_noise_i,sub_i,order_i,:) = interference_power;
                            
                            %                         SINR_optimal(est_i,speed_i,channel_noise_i,sub_i,order_i,:) = sinr_opt;
                            %mean_channel_matrix(speed_i,:,:,:) = mean_channel;
                            
                            %                         c = LTE_params.speed_of_light;
                            %                         f = LTE_params.carrier_freq;  % Frequency at which our system operates
                            %                         v = speed;  %speed at which we move
                            %                         w_d = 2*pi*v*f/c;   % maximum radian Doppler frequency
                            %                         T_s = LTE_params.SamplingTime*LTE_params.Nfft/LTE_params.Ntot;
                            %                         time_autocorrelation = besselj(0,w_d*(0:LTE_params.Ntot-1)*T_s);
                            
                            
                            % Code to generate the output filename
                            output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
                            filename_suffix = [];
                            
                            %save(fullfile('./results',[output_filename filename_suffix '.mat']));
                        end
                    end
                end
            end
        end
    end
end

figure(1)
hold on
grid on
plot(speed_vec*3.6,squeeze(tho(2,2,:,1,1,1)),'color','m','marker','d');
plot(speed_vec*3.6,squeeze(tho(2,2,:,1,1,5)),'color','g','marker','o');
plot(speed_vec*3.6,squeeze(tho(3,2,:,1,1,5)),'color','r','marker','s');
plot(speed_vec*3.6,squeeze(tho(4,2,:,1,1,5)),'color','b','marker','x');
xlabel('user speed [km/h]')
ylabel('throughput [Mbit/s]')

figure(2)
hold on
grid on
plot(speed_vec*3.6,squeeze(post_sinr(2,2,:,1,1,1)),'color','m','marker','d');
plot(speed_vec*3.6,squeeze(post_sinr(2,2,:,1,1,5)),'color','g','marker','o');
plot(speed_vec*3.6,squeeze(post_sinr(3,2,:,1,1,5)),'color','r','marker','s');
plot(speed_vec*3.6,squeeze(post_sinr(4,2,:,1,1,5)),'color','b','marker','x');
xlabel('user speed [km/h]')
ylabel('SINR [dB]')

figure(3)
hold on
grid on
plot(speed_vec*3.6,squeeze(tho(3,2,:,1,1,1)),'color','r','marker','s');
plot(speed_vec*3.6,squeeze(tho(4,2,:,1,1,1)),'color','r','marker','x');
plot(speed_vec*3.6,squeeze(tho(3,2,:,1,1,2)),'color','b','marker','s');
plot(speed_vec*3.6,squeeze(tho(4,2,:,1,1,2)),'color','b','marker','x');
plot(speed_vec*3.6,squeeze(tho(3,2,:,1,1,3)),'color','b','marker','s');
plot(speed_vec*3.6,squeeze(tho(4,2,:,1,1,3)),'color','b','marker','x');
xlabel('user speed [km/h]')
ylabel('throughput [Mbit/s]')

%save(fullfile('./results',[output_filename filename_suffix '.mat']));
% shutdown(10)