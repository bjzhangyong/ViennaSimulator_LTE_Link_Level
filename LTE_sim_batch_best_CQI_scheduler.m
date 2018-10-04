% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear
clear global
close all
% clc

%% DEBUG level
global DEBUG_LEVEL;
DEBUG_LEVEL = 4; % 1-5

Simulation_type = 'SUSISO_best_cqi';     %'SUSISO'
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
N_subframes = 15;
simulation_type = 'normal'; % 'parallel' or 'normal'
nr_re_tx_vect = [ 0 ];
cqi_i = 1;
for nr_re_tx = nr_re_tx_vect
    
    switch nr_re_tx
        case 0
            SNR_vec = 0:4:40;
        otherwise
            SNR_vec = -6:4:20;
    end
    
    % Load base parameter file
    LTE_load_parameters;
    LTE_params.scheduler.PMI = 0;
    LTE_params.max_HARQ_retransmissions = nr_re_tx;
    
    LTE_load_parameters_dependent;
    LTE_load_parameters_generate_elements;
    LTE_check_parameters;
    
    switch simulation_type
        case 'parallel'
            LTE_params.simulation_type = 'parallel';
        case 'normal'
            LTE_params.simulation_type = 'normal';
    end
    
    % The main simulation file. The rest of the code in this script is configuration/plotting
    LTE_sim_main;
    
    % Code to generate the output filename
    output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
    filename_suffix = [];
    
    save(fullfile('./results',[output_filename filename_suffix '.mat']));
    
end
