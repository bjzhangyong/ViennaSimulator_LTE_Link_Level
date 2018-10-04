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
SNR_end = [-4.1 -2.1 -0.1 1.9 3.9 5.9 5.9 9.9 11.9 13.9 15.9 17.4 18.9 20.65 22.4];
SNR_window = 20;
SNR_stepsize = 0.2;
SNR_begin = SNR_end-SNR_window;
N_subframes = 5000;

% Reorder simulations so that I get useful results faster
% cqi_list = 1:15;
cqi_list = 1:15;%[15 1 7 2 14 6 8 12 3 13 4 9 5 10 11];

Simulation_type = 'SUSISO_BLER_curves_batch';

%% Actual simulations
for cqi_i = cqi_list
    SNR_vec = SNR_begin(cqi_i):SNR_stepsize:SNR_end(cqi_i);
    LTE_load_parameters;  % BLER curves simulation
    
    % Change some parameters
    switch cqi_i
        case {1, 2, 3}
            LTE_params.Bandwidth = 20e6;
        case 4
            LTE_params.Bandwidth = 10e6;
        case {5, 6, 7}
            LTE_params.Bandwidth = 5e6;
        case {8, 9}
            LTE_params.Bandwidth = 3e6;
        otherwise
            LTE_params.Bandwidth = 1.4e6;
    end
    
    LTE_params.max_HARQ_retransmissions = 3;
    LTE_params.simulation_type = 'parallel';
    LTE_params.ChanMod_config.type = 'AWGN';
    LTE_params.random_noise_seeding = false;
    LTE_params.scheduler.type = 'round robin';
    LTE_params.scheduler.assignment = 'static';
    LTE_params.scheduler.cqi  = 'set';
    LTE_params.scheduler.PMI  = 0;
    
    % Redo config params
    LTE_load_parameters_dependent;
    LTE_load_parameters_generate_elements;
    LTE_check_parameters;
    
    LTE_sim_main;
    
    % Code to generate the output filename
    output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
    filename_suffix = [];
    
    save(fullfile('./results',[output_filename filename_suffix '.mat']));
end
