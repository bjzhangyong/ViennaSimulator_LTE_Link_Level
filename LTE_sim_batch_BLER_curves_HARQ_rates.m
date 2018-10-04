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
% 4-QAM (4 curves)
SNR_end(55) = 2;
SNR_end(56) = 3;
SNR_end(57) = 4;
SNR_end(58) = 5;
% 16-QAM (6 curves)
SNR_end(59) = 6;
SNR_end(60) = 6.5;
SNR_end(61) = 6.5;
SNR_end(62) = 7;
SNR_end(63) = 8;
SNR_end(64) = 10.5;
% 64-QAM (16 curves)
SNR_end(65) = 10.5;
SNR_end(66) = 10.5;
SNR_end(67) = 10.5;
SNR_end(68) = 10.5;
SNR_end(69) = 11;
SNR_end(70) = 12;
SNR_end(71) = 12;
SNR_end(72) = 12;
SNR_end(73) = 12.5;
SNR_end(74) = 13;
SNR_end(75) = 14;
SNR_end(76) = 14;
SNR_end(77) = 16;
SNR_end(78) = 17;
SNR_end(79) = 20;
SNR_end(80) = 22;

SNR_window = 5;
SNR_begin = SNR_end-SNR_window;
SNR_stepsize = 0.1;

% Reorder simulations so that I get useful results faster
% cqi_list = 1:15;
cqi_list = 55;

Simulation_type = 'SUSISO_BLER_curves_batch';

%% Plotting
first_figure = figure;
first_axes = axes;
set(first_axes,'YScale','log');
title('BLER, SISO AWGN, 5000 subframes');
ylabel('BLER');
xlabel('SNR [dB]');
ylim([1e-3 1]);
hold all;
grid on;

second_figure = figure;
second_axes =  axes;
title(second_axes,sprintf('throughput, SISO AWGN, 5000 subframes'));
ylabel(second_axes,'throughput [Mbps]');
xlabel(second_axes,'SNR [dB]');
hold all;
grid on;
                                    
%% Actual simulations
for cqi_i = cqi_list
    N_subframes = 5000;
    SNR_vec = SNR_begin(cqi_i):SNR_stepsize:SNR_end(cqi_i);
    LTE_load_parameters;  % BLER curves simulation
    
    % Change some parameters
    cqi_i_minus_54 = cqi_i-54;
    switch cqi_i_minus_54
        case {1}
            LTE_params.Bandwidth = 10e6;
        case {2,3,4,5,6,7}
            LTE_params.Bandwidth = 5e6;
        case {8,9,10,11,12,13,14,15,16}
            LTE_params.Bandwidth = 3e6;
        otherwise
            LTE_params.Bandwidth = 1.4e6;
    end
    LTE_params.max_HARQ_retransmissions = 0;
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
    
    % Plot
    name = sprintf('CQI %02d',cqi_i);
    plot(first_axes,simulation_results.SNR_vector,simulation_results.cell_specific.BLER_overall,'DisplayName',name);
    plot(second_axes,simulation_results.SNR_vector,squeeze(sum(sum(simulation_results.cell_specific.throughput_coded,1),3))/(N_subframes*1e-3)/1e6,'DisplayName',name);
end

legend(first_axes,'show','Location','SouthEastOutside');
legend(second_axes,'show','Location','SouthEastOutside');
