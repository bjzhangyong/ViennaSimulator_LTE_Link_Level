% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear
clear global
close all
clc

%% DEBUG level
global DEBUG_LEVEL;
DEBUG_LEVEL = 5; % Now set to highest level.

%% SNR setting
SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
SNR_stepsize = 0.1;
SNR_window = 3;

cqi_list = 1:15;

Simulation_type = 'SUSISO_BLER_curves_batch';     %'SUSISO'
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
for cqi_i = cqi_list
    N_subframes = 5000;
    SNR_vec = SNR_30percent(cqi_i)-SNR_window*2.5:SNR_stepsize:SNR_30percent(cqi_i)+SNR_window;
    
    LTE_load_parameters;  % Single User Single Input Single Output
    
    LTE_sim_main;
    
    % Code to generate the output filename
    output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
    filename_suffix = [];
    
    save(fullfile('./results',[output_filename filename_suffix '.mat']));
end
