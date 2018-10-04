classdef iaSimulationConfig < simulation_config.simulatorConfig
    methods (Static)
        function LTE_params = apply_parameters(LTE_params,N_Ue,N_Bs,channel_estimation_method,tx_mode,N_rx,N_tx,receiver,user_speed,filtering,channel_type,scheduler_type,scheduler_assignment,connection_table,Power_diff,IA_type,IA_streams,IA_thresh,IA_max_iterations,IA_freq_granularity,IA_time_granularity,IA_sigma_H2_E2_ratio)
            LTE_params.nUE = N_Ue;     % number of user equipments to simulate
            LTE_params.nBS = N_Bs;     % number of base stations to simulate (hard-coded to 1)
            LTE_params.Bandwidth = 1.4e6;            % in Hz, allowed values: 1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 15 MHz, 20MHz => number of resource blocks 6, 15, 25, 50, 75, 100
            LTE_params.introduce_timing_offset =  false;
            LTE_params.introduce_frequency_offset =  false;
            %% Define some User parameters (identical settings).
            LTE_params.UE_config.channel_estimation_method = channel_estimation_method;      %'PERFECT','LS','MMSE'
            LTE_params.UE_config.mode = tx_mode;                     % DEFINED IN STANDARD 3GPP TS 36.213-820 Section 7.1, page 12
            % 1: Single Antenna, 2: Transmit Diversity, 3: Open Loop Spatial Multiplexing
            % 4: Closed Loop SM, 5: Multiuser MIMO, 6: Interference Alignment
            LTE_params.UE_config.nRX = N_rx;                      % number of receive antennas at UE
            LTE_params.UE_config.receiver = receiver; % 'SSD','ZF'
            %         fd = 0; % Doppler frequency
            %         LTE_params.UE_config.user_speed = fd/LTE_params.carrier_freq*LTE_params.speed_of_light;  % [km/h]
            LTE_params.UE_config.user_speed = user_speed;
            LTE_params.UE_config.timing_offset = 23;   % timing offset in number of time samples
            LTE_params.UE_config.timing_sync_method = 'perfect';%'estimated';   % 'perfect','none','estimated','estimated2'
            LTE_params.UE_config.carrier_freq_offset = pi;   % carrier frequency offset normalized to subcarrier spacing
            LTE_params.UE_config.freq_sync_method = 'perfect';
            LTE_params.UE_config.rfo_correct_method = 'subframe'; % 'none','subframe'
            %         LTE_params.UE_config.user_speed = 300/3.6;    %[m/s]
            %         UE_params.SINR_averaging.averager = 'EESM';
            %% Define BS parameters (identical settings).
            LTE_params.BS_config.nTx = N_tx;
            %% Define ChanMod parameters - now it is only possible to have same channel parameters for BS and UE
            LTE_params.ChanMod_config.filtering = filtering;  %'BlockFading','FastFading'
            LTE_params.ChanMod_config.type = channel_type; % 'PedA', 'PedB', 'PedBcorr', 'AWGN', 'flat Rayleigh','VehA','VehB','TU','RA','HT','winner_II'
            %% Scheduler settings
            LTE_params.scheduler.type = scheduler_type;
            %        LTE_params.scheduler.type = 'constrained scheduler';
            % Available options are:
            %   - 'round robin': Will serve equally all of the available users
            %   - 'best cqi'   : Will serve only users that maximize the CQI for specific RB
            %   - 'fixed'
            
            LTE_params.scheduler.assignment = scheduler_assignment;
            %         LTE_params.scheduler.assignment = 'dynamic';
            % Available options are:
            %   - For 'round robin': 'static' of 'dynamic': whether the scheduler will statically
            %     assign or dynamically assign CQIs and other params. Currently only 'static' is implemented
            %   - For 'best cqi': 'dynamic': the scheduler will dynamically assign CQIs and other params.
            %   - For 'fixed': a vector stating how many RBs will each user get.
            
            % Parameters for the static scheduler
            LTE_params.scheduler.cqi  = 'set';
            LTE_params.scheduler.PMI  = 0;              % corresponds CI for closed loop SM
            
            LTE_params.connection_table = connection_table;
            LTE_params.Power_diff = Power_diff;
            LTE_params.Delay_diff = 0;  % 0 for synchronous, integer num of time samples
            
            LTE_params.IA_type = IA_type;
            LTE_params.IA_streams = IA_streams;
            LTE_params.IA_thresh = IA_thresh;
            LTE_params.IA_max_iterations = IA_max_iterations;
            LTE_params.IA_freq_granularity = IA_freq_granularity;
            LTE_params.IA_time_granularity = IA_time_granularity;
            LTE_params.IA_sigma_H2_E2_ratio = IA_sigma_H2_E2_ratio;
        end
    end
end

