% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear
clear global
close all
clc

%% DEBUG level
global DEBUG_LEVEL;
DEBUG_LEVEL = 4; % 1-5

%% SNR setting
SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
SNR_stepsize = 1;
SNR_window = 0.25;

speed_vec = [0:25:350]/3.6;
speed_vec(1) = 1/3.6;
tho = nan(size(speed_vec));

Simulation_type = 'winner_model_example';     %'SUSISO'
                                    %'MUSISO'
                                    %'SUMIMO'
                                    %'MUMIMO'
                                    %'SUSISO_quick_test'
                                    %'SUSISO_BLER_curves_batch'
                                    %'SUSISO_best_cqi'
                                    %'SUMIMO_quick_test'
                                    %'winner_model_example'
                                    %'wsa_2010_michal'
                                    
                                    
%% Actual simulations
for speed_i = 1:length(speed_vec)
    speed = speed_vec(speed_i)
    for cqi_i = 9
        N_subframes = 1000;
        SNR_vec = 20;
        
        LTE_load_parameters;  % Multi User Multiple Input Multiple Output using Winner II + Channel Model
        LTE_params.simulation_type = 'normal';

        % LTE_load_parameters_MUMIMO_winner_example;

        LTE_sim_main;
        tho(speed_i) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
        % Code to generate the output filename
        output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
        filename_suffix = [];

    %     save(fullfile('./results',[output_filename filename_suffix '.mat']));
    end
end
plot(speed_vec*3.6,tho);
xlabel('velocity [km/h]');
ylabel('Throughput [Mbit/s]');
save(fullfile('./results',[output_filename filename_suffix '.mat']));