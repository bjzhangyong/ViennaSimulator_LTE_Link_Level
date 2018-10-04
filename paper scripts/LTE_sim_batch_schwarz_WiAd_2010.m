% Batch file to reproduce the Figures shown in 
% "Calculation of the Spatial Preprocessing and Link Adaption Feedback for
% 3GPP UMTS/LTE" S. Schwarz, C. Mehlführer and M. Rupp
% NOTE: the optimal choice curves are not reproduced but loaded as the
% temporal effort for simulating this is very large. If you want to
% simulate them contact sschwarz@nt.tuwien.ac.at
% NOTE 2: Figure 8 cannot be reproduced because this feedback method is not
% used anymore in this simulator due to its poor performance
% If you experience any problems running the simulations contact sschwarz@nt.tuwien.ac.at



clear all
clear global
close all
clc
cd ..
%% DEBUG level
global DEBUG_LEVEL LTE_params;
DEBUG_LEVEL = 4;

Figure_choice = 'Fig1_2'; % change this to the Figure X from the paper you want to reproduce: FigX

%% Actual simulations
cqi_i = 4;
N_subframes = 5000;  
switch Figure_choice
    case 'Fig1_2' % Figures 1 and 2
        SNR_vec = -10:2:30;
        LTE_params.UE_config.nRX = 1;                      % number of receive antennas at UE
        LTE_params.BS_config.nTx = 2;
        LTE_params.UE_config.channel_estimation_method = 'PERFECT';      %'PERFECT','LS','MMSE'
        num_iter = 1;
        LTE_params.feedback.ignore_channel_estimation = false;
        LTE_params.feedback.channel_averaging = true;
    case 'Fig3' % Figure 3: just the optimal choice curves and the ones with estimated channel knowledge will be drawn. The curve with perfect channel knowledge can be obtained by choosing 'Fig1_2' 
        SNR_vec = -10:2:30;
        LTE_params.UE_config.nRX = 1;                      % number of receive antennas at UE
        LTE_params.BS_config.nTx = 2;
        LTE_params.UE_config.channel_estimation_method = 'LS';      %'PERFECT','LS','MMSE'
        num_iter = 2;
        LTE_params.feedback.channel_averaging = true;
    case 'Fig4'
        SNR_vec = -10:2:30;
        LTE_params.UE_config.nRX = 1;                      % number of receive antennas at UE
        LTE_params.BS_config.nTx = 2;
        LTE_params.UE_config.channel_estimation_method = 'MMSE';      %'PERFECT','LS','MMSE'
        num_iter = 2;
        LTE_params.feedback.channel_averaging = true;
    case 'Fig5_6'
        SNR_vec = -10:2:40;
        LTE_params.UE_config.nRX = 2;                      % number of receive antennas at UE
        LTE_params.BS_config.nTx = 2;
        LTE_params.UE_config.channel_estimation_method = 'PERFECT';      %'PERFECT','LS','MMSE'
        num_iter = 4;
        LTE_params.feedback.ignore_channel_estimation = false;
        LTE_params.feedback.channel_averaging = true;
    case 'Fig7'
        SNR_vec = -10:2:30;
        LTE_params.UE_config.nRX = 2;                      % number of receive antennas at UE
        LTE_params.BS_config.nTx = 4;
        LTE_params.UE_config.channel_estimation_method = 'PERFECT';      %'PERFECT','LS','MMSE'
        num_iter = 3;
        LTE_params.feedback.ignore_channel_estimation = false;
        LTE_params.feedback.channel_averaging = true;
end

for iterations = 1:num_iter
    if strcmp(Figure_choice,'Fig3') || strcmp(Figure_choice,'Fig4')
    if iterations == 1 
        LTE_params.feedback.ignore_channel_estimation = false;
    else
        LTE_params.feedback.ignore_channel_estimation = true;
    end
    end
%% Load necessary parameters
LTE_params.nUE = 1;     % number of user equipments to simulate
LTE_params.nBS = 1;     % number of base stations to simulate (hard-coded to 1)
LTE_params.Bandwidth = 1.4e6;            % in Hz, allowed values: 1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 15 MHz, 20MHz => number of resource blocks 6, 15, 25, 50, 75, 100
LTE_params.UE_config.mode = 4;                     % DEFINED IN STANDARD 3GPP TS 36.213-820 Section 7.1, page 12
LTE_params.UE_config.receiver = 'ZF'; % 'SSD','ZF'
LTE_params.UE_config.user_speed = 0/3.6;    % [km/h]
LTE_params.UE_config.freq_sync_method = 'perfect';
LTE_params.UE_config.rfo_correct_method = 'subframe'; % 'none','subframe'
LTE_params.ChanMod_config.filtering = 'BlockFading';  %'BlockFading','FastFading'
LTE_params.ChanMod_config.type = 'VehA'; % 'PedA', 'PedB', 'PedBcorr', 'AWGN', 'flat Rayleigh','VehA','VehB','TU','RA','HT','winner_II'
LTE_params.scheduler.type = 'fixed';
LTE_params.scheduler.assignment = 'semi static';
LTE_params.scheduler.cqi  = 'set';
LTE_params.scheduler.PMI  = 0;              % corresponds CI for closed loop SM
LTE_params.uplink_delay = 0; % Delay the uplink channel will introduce (in TTIs). >=1 TTI
LTE_params.show_plots = false;
LTE_params.plot_confidence_intervals = false;
LTE_params.trace_subcarrier_SNR = false;
LTE_params.N_seed_reset = 1;      % resets the random number generator seeds to a new value after LTE_params.N_seed_reset subframes
LTE_params.carrier_freq = 2.1e9;    %carrier frequency [Hz]
LTE_params.speed_of_light = 299792458; %[m/s]
LTE_params.HARQ_processes = 8;           % Number of HARQ processes
LTE_params.max_HARQ_retransmissions = 0; % max num of HARQ retransmissions, NOT including the original tx. 0, 1, 2 or 3
LTE_params.SubcarrierSpacing = 15e3; % in Hz, 15 kHz, also 7.5 kHz subcarrier spacing possible, just for MBSFN-based multicast/broadcast transmissions
LTE_params.CyclicPrefix = 'normal';  % 'normal' or 'extended' for MBSFN-based multicast/broadcast transmissions
LTE_params.simulation_type = 'parallel'; % 'parallel' or 'normal' to parallelize the SNR loop in LTE_sim_main_parallel.m for 4 cores or not as in LTE_sim_main.m
LTE_params.simulate_with_all_zero_sequences = false; % true if you want that the transmitted data is an all-zero sequence (useful for interleaver testing)
LTE_params.introduce_frequency_offset =  false;
LTE_params.random_noise_seeding = false; % Whether the seed for the random number generator that generates the noise is set
LTE_params.noise_seed = 0;               % Only used if the upper variable is set to 'true'
LTE_params.channel_matrix_source     = 'generated';  % 'generated' to generate every time, 'trace' to load it from a trace
LTE_params.store_channel_trace       = false;         % if mode is 'generated', the channel trace will be saved at the end of the simulation
LTE_params.channel_matrix_tracefile  = 'auto';       % filename where the trace is stored, only if 'trace' mode is used
LTE_params.CQI_mapping.coeffs = [0.5223 4.6176];    % these values hold for BLER 0.1 at first transmission
LTE_params.CQI_mapping.table = [-500;-6.934;-5.147;-3.18;-1.254;0.761;2.70;4.697;6.528;8.576;10.37;12.3;14.18;15.89;17.82;19.83;21]; % BLER 0.1 at first transmission
LTE_params.random_data_seeding = false;  % Whether the seed for the random number generator that generates the transmitted databits is set
LTE_params.data_seed = 10;   % Only used if the upper variable is set to 'true'
LTE_params.random_channel_param_seeding = false;  % Whether the seed for the random number generator that generates the channel parameters is set
LTE_params.channel_param_seed = 200;     % Only used if the upper variable is set to 'true'
LTE_params.UE_config.LLR_clipping = 100;
LTE_params.UE_config.turbo_iterations = 8;                       % Number of iterations of the turbo decoder
LTE_params.UE_config.N_soft = 1000000*LTE_params.HARQ_processes; % Defines the total number of soft channel bits available for HARQ processing (TS 36.306 4.2.1.3). Set to a high enough value.
LTE_params.UE_config.channel_interpolation_method = 'linear';    %'linear','cubic','spline','sinc_freq','sinc_time','T-F'
LTE_params.UE_config.autocorrelation_matrix_type 	= 'ideal';   % 'ideal','estimated'% type of autocorrelation amtrix('ideal','estimated')
LTE_params.UE_config.realization_num = 0;          % Number of realizations of channel, used for averaging fo channel autocorrelation matrix
LTE_params.UE_config.realization_num_total = 20;   % First xy number of channel realizations are used just for estimation of autocorrelation matrix
LTE_params.UE_config.CDD =0;                       % Cyclic Delay Diversity
LTE_params.UE_config.PMI_fb_granularity = 6;
LTE_params.UE_config.CQI_fb_granularity = 6;
LTE_params.UE_config.PMI_fb = true;             % PMI feedback activated (just used in connection with LTE_params.UE_config.mode = 4)
LTE_params.UE_config.RIandPMI_fb = true;           % RI and PMI feedback activated (just used in connection with LTE_params.UE_config.mode = 4)
LTE_params.UE_config.timing_offset = 23;   % timing offset in number of time samples
LTE_params.UE_config.timing_sync_method = 'perfect';% 'perfect','none', 'autocorrelation'
LTE_params.Pilots_power_offset = 0;
if strcmp('Fig5_6',Figure_choice) && (iterations == 3 || iterations == 4)
    LTE_params.UE_config.RIandPMI_fb = false;
end
if strcmp('Fig7',Figure_choice) && (iterations == 2 || iterations == 3)
    LTE_params.UE_config.RIandPMI_fb = false;
end
LTE_params.UE_config.CQI_fb = true;                 % CQI feedback activated 
LTE_params.UE_config.predict = false;               % channel prediction activated (used in the feedback calculation)
LTE_params.UE_config.carrier_freq_offset = pi;   % carrier frequency offset normalized to subcarrier spacing
LTE_params.UE_config.perfect_freq_sync = true;
LTE_params.UE_config.rfo_correct_method = 'none'; % 'none','subframe'
% LTE_params.UE_config.SINR_averaging.averager = 'EESM';    % SINR averaging method used for the feedback calculation (EESM or MIESM)
LTE_params.UE_config.SINR_averaging.averager = 'MIESM';    % SINR averaging method used for the feedback calculation (EESM or MIESM)
LTE_params.UE_config.SINR_averaging.EESMbetas = [5,5.01,5.01,0.84,1.67,1.61,1.64,3.87,5.06,6.4,12.59,17.59,23.33,29.45,33.05,35.41];  % weigthed EESM beta values
LTE_params.UE_config.SINR_averaging.MIESMbetas = [4,3.07,4.41,0.6,1.16,1.06,1.06,0.87,1.01,1.04,1.03,1.11,1.01,1.07,1,1.05];  % weighted MIESM beta values
LTE_params.UE_config.SINR_averaging.MCSs = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
LTE_params.scheduler.av_window = 10; % Scheduler averaging window size in subframes (not supported by all schedulers)
LTE_params.scheduler.fairness = 0.8; % desired fairness for the variable fairness scheduler (only supported by the var fair scheduler)
LTE_params.ChanMod_config.interpolation_method = 'shift_to_nearest_neighbor'; % 'shift_to_nearest_neighbor' 'sinc_interpolation'
LTE_params.ChanMod_config.corr_coefRX = 0.3;
LTE_params.ChanMod_config.corr_coefTX = 0.3;
LTE_params.ChanMod_config.sin_num = 10; % Number of sin realizations
LTE_params.ChanMod_config.time_correlation = 'independent';   % 'correlated','independent'
LTE_params.eig_tresh_time = 0.1;
LTE_params.eig_tresh_freq = 0.1;
LTE_params.ChanMod_config.winner_settings.Scenario                 = 11;                               % 1=A1, 2=A2, 3=B1, 4=B2, 5=B3, 6=B4, 7=B5a, 8=B5c, 9=B5f, 10=C1,
LTE_params.ChanMod_config.winner_settings.PropagCondition          = 'NLOS';                           % [LOS,{NLOS}]
LTE_params.ChanMod_config.winner_settings.SampleDensity            = 2;                                % number of time samples per half wavelength [ {2} ]
LTE_params.ChanMod_config.winner_settings.UniformTimeSampling      = 'yes';                            % Use same time sampling grid for all links [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.FixedPdpUsed             = 'no';            	               % nonrandom path delays and powers [ yes | {no}]
LTE_params.ChanMod_config.winner_settings.FixedAnglesUsed          = 'no';                             % nonrandom AoD/AoAs [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.PolarisedArrays          = 'yes';                            % usage of dual polarised arrays [ {yes} | no ]
LTE_params.ChanMod_config.winner_settings.TimeEvolution            = 'no';                             % usage of time evolution  [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.PathLossModelUsed        = 'no';                             % usage of path loss model [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.ShadowingModelUsed       = 'no';                             % usage of shadow fading model [ yes | {no} ]
LTE_params.ChanMod_config.winner_settings.PathLossModel            = 'pathloss';      	               % path loss model function name [ {pathloss} ]
LTE_params.ChanMod_config.winner_settings.PathLossOption           = 'CR_light';                       % ['{CR_light}' or 'CR_heavy' or 'RR_light' or 'RR_heavy', CR = Corridor-Room, RR = Room-Room nlos}
LTE_params.ChanMod_config.winner_settings.RandomSeed               = [];                               % sets random seed [ {[empty]} ]
LTE_params.ChanMod_config.winner_settings.UseManualPropCondition   = 'yes';                            % whether to use manual propagation condition (los/nlos) setting or not. If not, the propagation condition is drawn from probabilities.  [ {yes} | no]
LTE_params.usePBCH = false;
LTE_params.usePDCCH = false;
LTE_params.trafficmodel.usetraffic_model = false;

LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
LTE_params.Simulation_type = 'not_defined';
LTE_params.introduce_timing_offset = false;
LTE_params.SNR_estimation = false;

LTE_load_parameters_dependent;
LTE_load_parameters_generate_elements;
LTE_check_parameters;


if strcmp('Fig5_6',Figure_choice)
    switch iterations
        case 1
            LTE_params.feedback.channel_averaging = true;
        case 2
            LTE_params.feedback.channel_averaging = false;
        case 3
            LTE_params.feedback.channel_averaging = true;
            LTE_params.scheduler.nLayers(uu) =  2;
            LTE_params.scheduler.nCodewords(uu) =  min(2,LTE_params.scheduler.nLayers(uu));
        case 4
            LTE_params.feedback.channel_averaging = true;
            LTE_params.scheduler.nLayers(uu) =  1;
            LTE_params.scheduler.nCodewords(uu) =  min(2,LTE_params.scheduler.nLayers(uu));
    end
end
if strcmp('Fig7',Figure_choice)
    switch iterations
        case 1
            LTE_params.feedback.channel_averaging = true;
        case 2
            LTE_params.feedback.channel_averaging = true;
            LTE_params.scheduler.nLayers(uu) =  2;
            LTE_params.scheduler.nCodewords(uu) =  min(2,LTE_params.scheduler.nLayers(uu));
        case 3
            LTE_params.feedback.channel_averaging = true;
            LTE_params.scheduler.nLayers(uu) =  1;
            LTE_params.scheduler.nCodewords(uu) =  min(2,LTE_params.scheduler.nLayers(uu));
    end
end

%% Call the main program
LTE_sim_main
output_filename = ['Figures_' Figure_choice '_' num2str(iterations)];
save(fullfile([output_filename '.mat']));
end

switch Figure_choice
    case 'Fig1_2'
        figure(1)
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'rd-','Linewidth',2.5,'Markersize',14)
        hold on
        grid on
        set(gca,'Fontsize',24);
        set(gca,'Linewidth',1.5);
        xlabel('SNR [dB]');
        ylabel('Throughput [Mbit/s]');
        
        figure(2)
        semilogy(SNR_vec,simulation_results.UE_specific.BLER_overall,'bx-','Linewidth',2.5,'Markersize',14)
        set(gca,'Fontsize',24);
        set(gca,'Linewidth',1.5);
        xlabel('SNR [dB]');
        ylabel('BLER');
        
        load('examples/2x1_VehA_opt_choice.mat')
        figure(1)
        plot(SNR_vec,through_tot,'b+-','Linewidth',2.5,'Markersize',14)
        legend('MIESM feedback','Simulated optimal CQI, PMI');
   case 'Fig3'
        figure(1)
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'ms-','Linewidth',2.5,'Markersize',14)
        hold on
        grid on
        set(gca,'Fontsize',24);
        set(gca,'Linewidth',1.5);
        xlabel('SNR [dB]');
        ylabel('Throughput [Mbit/s]');
        load('Figures_Fig3_1.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'k+-','Linewidth',2.5,'Markersize',14)
        load('examples/2x1_VehA_opt_choice.mat')
        plot(SNR_vec,through_tot,'bx-','Linewidth',2.5,'Markersize',14)
        load('examples/2x1_VehA_LSestimator_optimal.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'go-','Linewidth',2.5,'Markersize',14)
        legend('CQI, PMI feedback, channel estimation error','CQI, PMI feedback, estimated channel','Simulated optimal CQI, PMI, perfect channel','Simulated optimal CQI, PMI, estimated channel','Location','Northwest')
    case 'Fig4'
        figure(1)
        load('Figures_Fig4_2.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'ms-','Linewidth',2.5,'Markersize',14)
        hold on
        grid on
        set(gca,'Fontsize',24);
        set(gca,'Linewidth',1.5);
        xlabel('SNR [dB]');
        ylabel('Throughput [Mbit/s]');
        load('Figures_Fig4_1.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'k+-','Linewidth',2.5,'Markersize',14)
        load('examples/2x1_VehA_opt_choice.mat')
        plot(SNR_vec,through_tot,'bx-','Linewidth',2.5,'Markersize',14)
        load('examples/VehA_2x1_optimal_MMSE.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'go-','Linewidth',2.5,'Markersize',14)
        legend('CQI, PMI feedback, channel estimation error','CQI, PMI feedback, estimated channel','Simulated optimal CQI, PMI, perfect channel','Simulated optimal CQI, PMI, estimated channel','Location','Northwest')
   case 'Fig5_6'
        figure(1)
        load('Figures_Fig5_6_1.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'g+-','Linewidth',2.5,'Markersize',14)
        hold on
        grid on
        set(gca,'Fontsize',24);
        set(gca,'Linewidth',1.5);
        xlabel('SNR [dB]');
        ylabel('Throughput [Mbit/s]');
        load('Figures_Fig5_6_2.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'ro-','Linewidth',2.5,'Markersize',14)
        load('Figures_Fig5_6_3.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'kd:','Linewidth',2.5,'Markersize',14)
        load('Figures_Fig5_6_4.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'ks:','Linewidth',2.5,'Markersize',14)
        load('examples/2x2_VehA_opt_PMI_CQI_RI.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'bo-','Linewidth',2.5,'Markersize',14)
        legend('CQI, PMI, RI feedback','CQI, PMI, RI feedback, no channel averaging','CQI, PMI feedback, RI = 2','CQI, PMI feedback, RI = 1','Optimal CQI, PMI, RI','Location','Northwest')
        
        figure(2)
        load('Figures_Fig5_6_1.mat');
        semilogy(SNR_vec,simulation_results.UE_specific.BLER_overall,'g+-','Linewidth',2.5,'Markersize',14)
        hold on
        grid on
        set(gca,'Fontsize',24);
        set(gca,'Linewidth',1.5);
        xlabel('SNR [dB]');
        ylabel('Throughput [Mbit/s]');
        load('Figures_Fig5_6_2.mat');
        semilogy(SNR_vec,simulation_results.UE_specific.BLER_overall,'ro-','Linewidth',2.5,'Markersize',14)
        load('Figures_Fig5_6_3.mat');
        semilogy(SNR_vec,simulation_results.UE_specific.BLER_overall,'kd:','Linewidth',2.5,'Markersize',14)
        load('Figures_Fig5_6_4.mat');
        semilogy(SNR_vec,simulation_results.UE_specific.BLER_overall,'ks:','Linewidth',2.5,'Markersize',14)
%         load('examples/2x2_VehA_opt_PMI_CQI_RI.mat');
%         semilogy(SNR_vec,simulation_results.UE_specific.BLER_overall,'bo-','Linewidth',2.5,'Markersize',14)
        legend('CQI, PMI, RI feedback','CQI, PMI, RI feedback, no channel averaging','CQI, PMI feedback, RI = 2','CQI, PMI feedback, RI = 1','Location','Northwest')
   case 'Fig7'
        figure(1)
        h = zeros(4,1);
        load('Figures_Fig7_1.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'g+-','Linewidth',2.5,'Markersize',14);
        hold on
        grid on
        set(gca,'Fontsize',24);
        set(gca,'Linewidth',1.5);
        xlabel('SNR [dB]');
        ylabel('Throughput [Mbit/s]');
        load('Figures_Fig7_2.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'ro-','Linewidth',2.5,'Markersize',14);
        load('Figures_Fig7_3.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'kd:','Linewidth',2.5,'Markersize',14);
        load('examples/4x2_VehA_opt_PMI_CQI_RI.mat');
        plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'bo-','Linewidth',2.5,'Markersize',14);
        legend('CQI, PMI, RI feedback','CQI, PMI feedback, RI = 2','CQI, PMI feedback, RI = 1','Optimal CQI, PMI, RI','Location','Northwest')
 end
