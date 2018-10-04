% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

%clear all
%clear global
close all
clc

%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;

%% SNR setting
% SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
% SNR_stepsize = 1;
% SNR_window = 0.25;
% power = [];
% noise = [];
Simulation_type = 'SUSISO_best_cqi';
SINR_data = load('./paper scripts/HARQ_validation_SINRs.mat');
the_SNR_vector = SINR_data.eNodeB_sector_SINRs_cut';

%% Actual simulations
SNR_vec = the_SNR_vector;

%SNR_vec = the_SNR_vector;

N_subframes = 500;
cqi_i = 0; %(dummy parameter)
LTE_load_parameters;  % BLER curves simulation

% Change some parameters
LTE_params.nUE = 1;
LTE_params.max_HARQ_retransmissions = 3;
LTE_params.simulation_type = 'parallel';
LTE_params.ChanMod_config.type = 'PedB';
LTE_params.UE_config.user_speed = 5/3.6;
LTE_params.random_noise_seeding = false;
LTE_params.scheduler.type = 'best cqi HARQ single-user';
LTE_params.scheduler.cqi  = 'dynamic';
LTE_params.scheduler.PMI  = 0;
LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);

% Redo config params
LTE_load_parameters_dependent;
LTE_load_parameters_generate_elements;
LTE_check_parameters;

LTE_sim_main;

% Code to generate the output filename
output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
filename_suffix = [];

save(fullfile('./results',[output_filename filename_suffix '.mat']));

