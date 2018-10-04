% Configuration file for Multi-User MIMO
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2008/08/11
% last changes: 2008/09/02  Colom Ikuno new variables for HARQ
%               2008/09/15  Bosanska    changed structure of BS and UE
%               2008/09/18  Bosanska    added definition of the channel model [nBS x nUE]struct ChanMod
%                                       added Initial state of random number generator
%               2008/10/21  Bosanska    added LTE_params.config.CQI_model -> linear regression for the CQI mapping
%               2008/11/05  Bosanska    extended structure ChanMod to design of RRC TX and RX filters
%               2008/11/07  Bosanska    removed extended structure ChanMod to design of RRC TX and RX filters
%                                       channel taps for PedA and PedB are given in [s], are shifted according to the
%                                       sampling freq LTE_params.Fs and rounded
%               2008/12/10  Simko       added type of channel estimation
%                                       method to UE -> UE.channel_estimation_method
%
% Contains code the following sources:
%   - Coded Modulation Library (CML) by Iterative Solutions:

global LTE_params;

LTE_params.Simulation_type = Simulation_type;

LTE_params.carrier_freq = 2.5e9;    %carrier frequency [Hz]
LTE_params.speed_of_light = 299792458; %[m/s]
switch LTE_params.Simulation_type
    case 'SUSISO'
        LTE_params = simulation_config.susisoSimulationConfig.apply_parameters(LTE_params);
    case 'MUSISO'
        LTE_params = simulation_config.musisoSimulationConfig.apply_parameters(LTE_params);
    case 'SUMIMO'
        LTE_params = simulation_config.sumimoSimulationConfig.apply_parameters(LTE_params);
    case 'MUMIMO'
        LTE_params = simulation_config.mumimoSimulationConfig.apply_parameters(LTE_params);
    case 'wsa_2010_michal'
        LTE_params = simulation_config.wsa2010michalSimulationConfig.apply_parameters(LTE_params,estimator,speed);
    case 'wsa_2010_schwarz'
        LTE_params = simulation_config.wsa2010schwarzSimulationConfig.apply_parameters(LTE_params,nRX,nTX,receiver,channel);
    case 'pimrc_2010_qwang'
        LTE_params = simulation_config.pimrc2010qwangSimulationConfig.apply_parameters(LTE_params,CFO,freq_sync_method,channel_type,mode,nRX,nTX);
    case 'winner_model_example'
        LTE_params = simulation_config.winnerModelExampleSimulationConfig.apply_parameters(LTE_params);
    case 'SUSISO_quick_test'
        LTE_params = simulation_config.sisoQuickTestSimulationConfig.apply_parameters(LTE_params);
    case 'SUSISO_BLER_curves_batch'
        LTE_params = simulation_config.susisoBlerCurvesBatchSimulationConfig.apply_parameters(LTE_params);
    case 'SUSISO_best_cqi'
        LTE_params = simulation_config.susisoBestCqiSimulationConfig.apply_parameters(LTE_params);
    case 'SUMIMO_quick_test'
        LTE_params = simulation_config.sumimoQuickTestSimulationConfig.apply_parameters(LTE_params);
    case 'TB'
        LTE_params = simulation_config.tbSimulationConfig.apply_parameters(LTE_params,Testbed_settings);
    case 'beta'
        LTE_params = simulation_config.betaSimulationConfig.apply_parameters(LTE_params,channel_type);
    case 'LTE_journal_paper'
        LTE_params = simulation_config.lteJournalPaperSimulationConfig.apply_parameters(LTE_params,N_Ue,N_Bs,tx_mode,N_rx,N_tx,channel_type,connection_table,Power_diff,channel_est,equalizer);
    case 'vtc_2011_spring_michal'
        LTE_params = simulation_config.vtc2011springMichalSimulationConfig.apply_parameters(LTE_params,estimator,speed);
    case 'IA'
        LTE_params = simulation_config.iaSimulationConfig.apply_parameters(LTE_params,N_Ue,N_Bs,channel_estimation_method,tx_mode,N_rx,N_tx,receiver,user_speed,filtering,channel_type,scheduler_type,scheduler_assignment,connection_table,Power_diff,IA_type,IA_streams,IA_thresh,IA_max_iterations,IA_freq_granularity,IA_time_granularity,IA_sigma_H2_E2_ratio);    
    case 'mobilkom_CLSM'
        LTE_params = simulation_config.mobilkom_CLSM.apply_parameters(LTE_params);
    otherwise
        error('not valid Simulation type')
end

% %% Example parameters that are written by the prior lines (SUSISO example)
% LTE_params.nUE = 1;     % number of user equipments to simulate
% LTE_params.nBS = 1;     % number of base stations to simulate (hard-coded to 1)
% LTE_params.Bandwidth = 1.4e6;            % in Hz, allowed values: 1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 15 MHz, 20MHz => number of resource blocks 6, 15, 25, 50, 75, 100
% LTE_params.introduce_frequency_offset =  false;
% %% Define some User parameters (identical settings).
% LTE_params.UE_config.channel_estimation_method = 'LS';      %'PERFECT','LS','MMSE'
% LTE_params.UE_config.mode = 1;                     % DEFINED IN STANDARD 3GPP TS 36.213-820 Section 7.1, page 12
% % 1: Single Antenna, 2: Transmit Diversity, 3: Open Loop Spatial Multiplexing
% % 4: Closed Loop SM, 5:
% % Multiuser MIMO
% LTE_params.UE_config.nRX = 1;                      % number of receive antennas at UE
% LTE_params.UE_config.receiver = 'ZF'; % 'SSD','ZF'
% LTE_params.UE_config.user_speed = 10/3.6;    %[m/s]
% LTE_params.UE_config.carrier_freq_offset = pi;   % carrier frequency offset normalized to subcarrier spacing
% LTE_params.UE_config.freq_sync_method = 'perfect';
% LTE_params.UE_config.rfo_correct_method = 'subframe'; % 'none','subframe'
% %% Define BS parameters (identical settings).
% LTE_params.BS_config.nTx = 1;
% %% Define ChanMod parameters - now it is only possible to have same channel parameters for BS and UE
% LTE_params.ChanMod_config.filtering = 'BlockFading';  %'BlockFading','FastFading'
% LTE_params.ChanMod_config.type = 'PedB'; % 'PedA', 'PedB', 'PedBcorr', 'AWGN', 'flat Rayleigh','VehA','VehB','TU','RA','HT','winner_II'
% %% Scheduler settings
% LTE_params.scheduler.type = 'round robin';
% % Available options are:
% %   - 'round robin': Will serve equally all of the available users
% %   - 'best cqi'   : Will serve only users that maximize the CQI for specific RB
% %   - 'fixed'
% 
% LTE_params.scheduler.assignment = 'static';
% % Available options are:
% %   - For 'round robin': 'static' of 'dynamic': whether the scheduler will statically
% %     assign or dynamically assign CQIs and other params. Currently only 'static' is implemented
% %   - For 'best cqi': 'dynamic': the scheduler will dynamically assign CQIs and other params.
% %   - For 'fixed': a vector stating how many RBs will each user get.
% 
% % Parameters for the static scheduler
% LTE_params.scheduler.cqi  = 'set';
% LTE_params.scheduler.PMI  = 0;              % corresponds CI for closed loop SM
% %% End of the SUSISO example

if ~exist('connection_table','var')
    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
else
    LTE_params.connection_table = connection_table;
end

if ~exist('delay_table','var')
    LTE_params.delay_table = zeros(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
else
    LTE_params.delay_table = zeros(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);  % asynchronous network case, change here.
end


LTE_params.SNR_estimation = false; %true-real SNR estimation, false-perfect SNR knowledge
LTE_params.number_of_zeros_for_SNR_estimation = 256;
if strcmp('TB',LTE_params.Simulation_type)
    LTE_params.number_of_zeros_for_SNR_estimation = Testbed_settings.Nr_of_nulls;
end

LTE_params.uplink_delay = 1; % Delay the uplink channel will introduce (in TTIs). >=1 TTI
LTE_params.show_plots = false;
LTE_params.plot_confidence_intervals = true;
LTE_params.trace_subcarrier_SNR = false;
LTE_params.use_seed_reset = false; %wheter the channel seed should be reseted or not
LTE_params.N_seed_reset = 1;      % resets the random number generator seeds to a new value after LTE_params.N_seed_reset subframes
if ~exist('pilot_offset','var')
    LTE_params.Pilots_power_offset = 0; % [dB] offset of the pilot energy compared to teh data while keeping the total transmit energy the same
else
    LTE_params.Pilots_power_offset = pilot_offset;
end
%PRIMITIVE PARAMETERS
%DEFINED IN STANDARD 3GPP TS 36.104 V8.1.0 (2008-03), page 10
%Table 5.1-1 Transmission bandwidth configuration
LTE_params.HARQ_processes = 8;           % Number of HARQ processes
LTE_params.max_HARQ_retransmissions = 0; % max num of HARQ retransmissions, NOT including the original tx. 0, 1, 2 or 3

LTE_params.SubcarrierSpacing = 15e3; % in Hz, 15 kHz, also 7.5 kHz subcarrier spacing possible, just for MBSFN-based multicast/broadcast transmissions
LTE_params.CyclicPrefix = 'normal';  % 'normal' or 'extended' for MBSFN-based multicast/broadcast transmissions

LTE_params.simulation_type = 'normal'; % 'parallel' or 'normal' to parallelize the SNR loop in LTE_sim_main_parallel.m for 4 cores or not as in LTE_sim_main.m
LTE_params.simulate_with_all_zero_sequences = false; % true if you want that the transmitted data is an all-zero sequence (useful for interleaver testing)

%% Configuration of the random noise generation
LTE_params.random_noise_seeding = false; % Whether the seed for the random number generator that generates the noise is set
LTE_params.noise_seed = 0;               % Only used if the upper variable is set to 'true'

%% Channel matrix source
LTE_params.channel_matrix_source     = 'generated';  % 'generated' to generate every time, 'trace' to load it from a trace
LTE_params.store_channel_trace       = false;         % if mode is 'generated', the channel trace will be saved at the end of the simulation
LTE_params.channel_matrix_tracefile  = 'auto';       % filename where the trace is stored, only if 'trace' mode is used

%% Define config parameters for specific algorithms
% CQI mapping
% LTE_params.CQI_mapping.coeffs = [0.5298 4.3002]; % linear regression for the SINR to CQI mapping
LTE_params.CQI_mapping.coeffs = [0.5223 4.6176];    % these values hold for BLER 0.1 at first transmission
LTE_params.CQI_mapping.table = [-500;-6.934;-5.147;-3.18;-1.254;0.761;2.70;4.697;6.528;8.576;10.37;12.3;14.18;15.89;17.82;19.83;21]; % BLER 0.1 at first transmission
% LTE_params.CQI_mapping.table = [-500;-7;-5.661;-3.595;-1.572;0.537;2.545;4.593;6.408;8.479;10.34;12.22;14.13;15.85;17.79;19.86;22]; % max throughput AWGN
% LTE_params.CQI_mapping.table = [-500;-6.936;-5.147;-3.121;-1.175;0.815;2.63;4.67;6.5;8.5;10.6;12.5;15.5;16.1;17.8;19.7;22]; % 4x2 taylored
% LTE_params.CQI_mapping.table = [-500;-6.936;-5.147;-3.121;-1.175;0.761;2.63;4.67;6.5;8.5;10.4;12;14.5;16.1;18;22;24]; % 4x4 taylored

%% Initial state of random number generator
LTE_params.random_data_seeding = true;  % Whether the seed for the random number generator that generates the transmitted databits is set
LTE_params.data_seed = 10;   % Only used if the upper variable is set to 'true'
LTE_params.random_channel_param_seeding = true;  % Whether the seed for the random number generator that generates the channel parameters is set
LTE_params.channel_param_seed = 175;     % Only used if the upper variable is set to 'true'

%% Define UE parameters (the same for all users, although it could be changed to have different users)
LTE_params.UE_config.LLR_clipping = 100;
LTE_params.UE_config.decoder_type = 'max-log-map';
LTE_params.UE_config.turbo_iterations = 8;                       % Number of iterations of the turbo decoder
LTE_params.UE_config.hard_demapping = false;
LTE_params.UE_config.N_soft = 3667200;                           % Defines the total number of soft channel bits available (TS 36.306, 4.2.1.3). Set to category 5 (Table 4.1-1)
LTE_params.UE_config.channel_interpolation_method = 'linear';    %'linear','cubic','spline','sinc_freq','sinc_time','T-F'
% interpolation for fast fading 'linear','cubic','v4'
LTE_params.UE_config.autocorrelation_matrix_type 	= 'ideal';   % 'ideal','estimated'% type of autocorrelation amtrix('ideal','estimated')
LTE_params.UE_config.realization_num = 0;          % Number of realizations of channel, used for averaging fo channel autocorrelation matrix
LTE_params.UE_config.realization_num_total = 20;   % First xy number of channel realizations are used just for estimation of autocorrelation matrix
LTE_params.UE_config.CDD =0;                       % Cyclic Delay Diversity
% 0... zero delay CDD (3GPP TS 36.211-820 Section 6.3.4.2.1, page 37)
% 1... small delay CDD (3GPP TS 36.211-820 Section 6.3.4.2.1, page 37) --> in the newest standard version this is not defined anymore
% 2... large delay CDD (3GPP TS 36.211-820 Section 6.3.4.2.2, page 38)

switch LTE_params.Bandwidth                         % choose coarsest PMI feedback granularity (or use the line below)
    case 1.4*10^6
        LTE_params.UE_config.PMI_fb_granularity = 6;
        LTE_params.UE_config.CQI_fb_granularity = 6;
    case 3*10^6
        LTE_params.UE_config.PMI_fb_granularity = 15;
        LTE_params.UE_config.CQI_fb_granularity = 15;
    case 5*10^6
        LTE_params.UE_config.PMI_fb_granularity = 25;
        LTE_params.UE_config.CQI_fb_granularity = 25;
    case 10*10^6
        LTE_params.UE_config.PMI_fb_granularity = 50;
        LTE_params.UE_config.CQI_fb_granularity = 50;
    case 15*10^6
        LTE_params.UE_config.PMI_fb_granularity = 75;
        LTE_params.UE_config.CQI_fb_granularity = 75;
    case 20*10^6
        LTE_params.UE_config.PMI_fb_granularity = 100;
        LTE_params.UE_config.CQI_fb_granularity = 100;
end
% LTE_params.UE_config.PMI_fb_granularity = 1;       % PMI feedback granularity in multiples of resource blocks (only one value per full bandwidth or per resource block is supported)
% LTE_params.UE_config.CQI_fb_granularity = 1;       % CQI feedback granularity in multiples of RBs (only one value per full bandwidth or per resource block is supported)
if strcmp(LTE_params.Simulation_type,'wsa_2010_schwarz')
    LTE_params.UE_config.PMI_fb = PMI_fb;
    LTE_params.UE_config.RIandPMI_fb = false;
    LTE_params.UE_config.CQI_fb = false;
    LTE_params.UE_config.PMI_fb_granularity = 6;
    LTE_params.UE_config.CQI_fb_granularity = 6;
else
    LTE_params.UE_config.PMI_fb = true;             % PMI feedback activated (just used in connection with LTE_params.UE_config.mode = 4)
    LTE_params.UE_config.RIandPMI_fb = true;           % RI and PMI feedback activated (just used in connection with LTE_params.UE_config.mode = 4)
    LTE_params.UE_config.CQI_fb = true;                 % CQI feedback activated
end
LTE_params.UE_config.predict = false;               % channel prediction activated (used in the feedback calculation)

% LTE_params.UE_config.SINR_averaging.averager = 'EESM';    % SINR averaging method used for the feedback calculation (EESM or MIESM)
LTE_params.UE_config.SINR_averaging.averager = 'MIESM';    % SINR averaging method used for the feedback calculation (EESM or MIESM)
% LTE_params.UE_config.SINR_averaging.EESMbetas = [4,3.96,2.77,0.93,1.38,1.44,1.56,3.62,4.73,6.09,11.59,16.35,21.45,26.5,30.75,33.5];  % EESM beta values
% LTE_params.UE_config.SINR_averaging.MIESMbetas = [4,3.85,2.85,0.66,1.04,0.98,1,0.82,0.95,1,0.99,1.02,0.94,1.03,1,1];  % MIESM beta values
LTE_params.UE_config.SINR_averaging.EESMbetas = [5,5.01,5.01,0.84,1.67,1.61,1.64,3.87,5.06,6.4,12.59,17.59,23.33,29.45,33.05,35.41];  % weigthed EESM beta values
LTE_params.UE_config.SINR_averaging.MIESMbetas = [4,3.07,4.41,0.6,1.16,1.06,1.06,0.87,1.01,1.04,1.03,1.11,1.01,1.07,1,1.05];  % weighted MIESM beta values

% LTE_params.UE_config.SINR_averaging.MIESMbetas = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];  % MIESM beta values
LTE_params.UE_config.SINR_averaging.MCSs = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];

%% Other parameters
% Parameters for the feedback calculation
LTE_params.feedback.ignore_channel_estimation = false; % ignores the channel estimation MSE for the feedback calculation if set
LTE_params.feedback.channel_averaging = true; % use channel averaging for feedback calculation

% Scheduler parameters
LTE_params.scheduler.av_window = 100; % Scheduler averaging window size in subframes (not supported by all schedulers) - 100 is recommended for the varfair algorithm
LTE_params.scheduler.fairness = 0.95; % desired fairness for the variable fairness scheduler (only supported by the 'var fair' scheduler)
LTE_params.scheduler.alpha = 0;   % alpha level for the 'alpha fair' scheduler
LTE_params.scheduler.weights = 1/LTE_params.nUE*ones(LTE_params.nUE,1); % weights for the weighted sum rate maximization in the 'weighted sum rate' scheduler
LTE_params.scheduler.stepsize = 0.5; % initial stepsize for the adaptation of the parameter alpha in the startup phase of the 'var fair' scheduler - 0.5 is recommended
LTE_params.scheduler.rate_constraints = Inf(LTE_params.nUE,1);  % average rate constraints applied in the utility maximizing scheduler
LTE_params.scheduler.rate_constraints = [Inf;Inf];
LTE_params.scheduler.MUMIMO = false; % whether MUMIMO is used with the constrained scheduler or SUMIMO

%% Define ChanMod parameters - now it is only possible to have same channel parameters for BS and UE
LTE_params.ChanMod_config.interpolation_method = 'shift_to_nearest_neighbor'; % 'shift_to_nearest_neighbor' 'sinc_interpolation'
LTE_params.ChanMod_config.corr_coefRX = 0.3;
LTE_params.ChanMod_config.corr_coefTX = 0.3;
LTE_params.ChanMod_config.sin_num = 10; % Number of sin realizations
LTE_params.ChanMod_config.time_correlation = 'independent';   % 'correlated','independent'
LTE_params.eig_tresh_time = 0.1;
LTE_params.eig_tresh_freq = 0.1;
% winner model settings
LTE_params.ChanMod_config.winner_settings.Scenario                 = 1; %10                            % 1=A1, 2=A2, 3=B1, 4=B2, 5=B3, 6=B4, 7=B5a, 8=B5c, 9=B5f, 10=C1,
% 11=C2, 12=C3, 13=C4, 14=D1 and 15=D2a
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

%% include other channels (PBCH, PDCCH)
LTE_params.usePBCH = false;
LTE_params.usePDCCH = false;

%% use traffic models for simulation (true) or fullbuffer assumption (false) (note: this is just used if the scheduler supports it!)
LTE_params.trafficmodel.usetraffic_model = false; % just supported by the constrained scheduler

%% Load and calculate other parameters which depend on the ones previously defined
LTE_load_parameters_dependent;
LTE_load_parameters_generate_elements;
LTE_check_parameters;

