% Basic batch simulation script to recreate the BLER plots in the R1-07196 RAN document
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
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
SNR_mins = [-7 -6 -5 -4 -3 -2 -1 0 2 3 4 4 5 6 6 7 7 8 10 11 12 12 13 14 15 16 17];
SNR_maxs = SNR_mins+4;
SNR_stepsize = 0.25;

% Simulation prepared for 4 cores
core = 1;
N_subframes = 5000;
simulation_type = 'parallel'; % 'parallel' or 'normal'

switch core
    case 1
        cqis = [101:103 127:-1:124];
    case 2
        cqis = [104:106 123:-1:120];
    case 3
        cqis = [107:109 119:-1:116];
    case 4
        cqis = [110:115];
end

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
for cqi_i = cqis
    SNR_vec = SNR_mins(cqi_i-100):SNR_stepsize:SNR_maxs(cqi_i-100);
    
    LTE_load_parameters;  % Single User Single Input Single Output
    
    switch simulation_type
        case 'parallel'
            LTE_params.simulation_type = 'parallel';
        case 'normal'
            LTE_params.simulation_type = 'normal';
    end
    
    % See comments in LTE_sim_main for using parfor
    LTE_sim_main;

    % Code to generate the output filename
    output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes)
    filename_suffix = [];

    save(fullfile('./results',[output_filename filename_suffix '.mat']));
    %close all;
end
