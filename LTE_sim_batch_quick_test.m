% Script that gives you an overview of how to use the LTE link level
% simulator for a simple comparison test of the several TX modes in LTE.
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

% Number of subframes to simulate
N_subframes = 5000;
simulation_type = 'parallel'; % 'parallel' or 'normal'
channel_type = 'PedB';

% Example configuration for CQI 7
cqis_to_simulate = 7;
tx_modes_to_simulate = [111 221 242 342]; % XYZ--> TX mode, nTX antennas, nRX antennas
nr_re_tx_vect = [ 0 3 ];

MIMO         = tx_modes_to_simulate >= 200;
tx_modes     = floor(tx_modes_to_simulate/100);
nTX_antennas = floor(mod(tx_modes_to_simulate,100)/10);
nRX_antennas = mod(mod(tx_modes_to_simulate,100),10);

% NOTE: not all combinatios are defined in this script

for nr_re_tx = nr_re_tx_vect
    
    switch nr_re_tx
        case 0
            SNR_vec = 0:30;
        otherwise
            SNR_vec = -6:20;
    end
    
    for cqi_i = cqis_to_simulate
        
        % Actual simulations
        for simulation_idx = 1:length(tx_modes_to_simulate)
            
            if MIMO(simulation_idx)
                Simulation_type = 'SUMIMO_quick_test';
            else
                Simulation_type = 'SUSISO_quick_test';
            end
            LTE_load_parameters;
            % Perform some changes
            LTE_params.UE_config.nRX  = nRX_antennas(simulation_idx); %Actually this is not necessary, but it makes it more clear
            LTE_params.BS_config.nTx  = nTX_antennas(simulation_idx); % Also not necessary, but imprves readability
            LTE_params.UE_config.mode = tx_modes(simulation_idx);
            LTE_params.show_plots = false;
            LTE_params.ChanMod_config.type = channel_type;
            LTE_params.trace_subcarrier_SNR = false; % This takes a lot of disk space, so if you don't need it, better to leave it off.
            LTE_params.max_HARQ_retransmissions = nr_re_tx;
            LTE_params.scheduler.cqi = 'set';
            LTE_params.scheduler.PMI  = 0;  
            % This file reads the configurations form before and
            % generates the needed structures for the program to run
            LTE_load_parameters_dependent;
            LTE_load_parameters_generate_elements
            % A file that does some consistency checks
            LTE_check_parameters;
            
            switch tx_modes(simulation_idx)
                case 1
                    name = 'SISO';
                case 2
                    name = sprintf('TxD %dx%d',nTX_antennas(simulation_idx),nRX_antennas(simulation_idx));
                case 3
                    name = sprintf('OLSPM %dx%d',nTX_antennas(simulation_idx),nRX_antennas(simulation_idx));
                case 4
                    name = sprintf('CLSPM %dx%d',nTX_antennas(simulation_idx),nRX_antennas(simulation_idx));
                otherwise
                    error('extend the tests!!!');
            end
            
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
    end
end
