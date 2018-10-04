classdef UE < handle
% Class that represents an LTE UE.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

properties
    demapping_method                % Demapping method to use
    LLR_clipping                    % LLR clipping
    turbo_iterations                % Number of iterations of the turbo decoder
    N_soft                          % Defines the total number of soft channel bits available for HARQ processing (TS 36.306 4.2.1.3)
    channel_estimation_method       % Channel estimation method
    channel_interpolation_method    % channel interpolation method
    channel_autocorrelation_matrix  % channel autocorrelation matrix
    autocorrelation_matrix_type 	% type of autocorrelation amtrix('ideal','estiamted')
    realization_num                 % number of realizations of channel, used for averaging fo channel autocorrelation matrix
    realization_num_total           % first xy number of channel realizations are used just for estimation of autocorrelation matrix
    user_speed                      % paramter for fast fading channel
    HARQ_rx_soft_buffer             % HARQ soft buffer. Stores codewords, not data bits!!!!
    CDD                             % Cyclic Delay Diversity
                                    % 0... zero delay CDD (3GPP TS 36.211-820 Section 6.3.4.2.1, page 37)
                                    % 1... small delay CDD (3GPP TS 36.211-820 Section 6.3.4.2.1, page 37)
                                    % 2... large delay CDD (3GPP TS 36.211-820 Section 6.3.4.2.2, page 38)
    mode                            % DEFINED IN STANDARD 3GPP TS 36.213-820 Section 7.1, page 12
                                    % 1: Single Antenna, 2: Transmit Diversity, 3: Open Loop Spatial Multiplexing
                                    % 4: Closed Loop SM, 5: Multiuser MIMO
    nRX                             % number of receive antennas at UE
    clock = 0;                      % So the UE is aware in which TTI he is (initialized to 0)
    timing_offset;                  % symbol timing offset
    timing_sync_method              % timing synchronization method
    carrier_freq_offset             % carrier frequency offset
    freq_sync_method                % carrier frequency offset correction method
    rfo_correct_method              % residual frequency offset compensation method
    receiver                        % receiver that the UE applys
    receiver_k                      % number of survivor nodes for k-best receiver
    channel_coef_rosa_zheng         % phi, theta, psi for Rosa Zhneg Channel Model
    previous_channels               % necessary for channel extrapolation
    PMI_fb_gran                     % PMI feedback granularity in multiples of resource blocks
    PMI_fb                          % wheter PMI feedback is activated
    RIandPMI_fb                     % wheter RI and PMI feedback is activated
    predict                         % wheter channel prediction is activated
    CQI_fb_gran                     % CQI feedback granularity 
    CQI_fb                          % wheter CQI feedback is activated
    SINR_averager                   % defines the SINR averaging method used for feedback calculation
    MSE                             % channel estimation MSE
    traffic_model                   % traffic model used for this UE
end

   methods
       % Class constructor, defining default values.
       % Get as input a struct defining the input variables
       function obj = UE(UE_params,Ntot,Nsub,nTX,fading)
           
           % Default values
           if ~isfield(UE_params,'turbo_iterations')
               UE_params.turbo_iterations = 8;
           end
           if ~isfield(UE_params,'rfo_correct_method')
               UE_params.rfo_correct_method = 'none';
           end
           if ~isfield(UE_params,'receiver_k')
               UE_params.receiver_k = 8;
           end

           % Assign values
           % obj.demapping_method             = UE_params.demapping_method;
           obj.LLR_clipping                 = UE_params.LLR_clipping;
           obj.turbo_iterations             = UE_params.turbo_iterations;
           obj.N_soft                       = UE_params.N_soft;
           obj.channel_estimation_method    = UE_params.channel_estimation_method;
           obj.channel_interpolation_method = UE_params.channel_interpolation_method;
           obj.autocorrelation_matrix_type 	= UE_params.autocorrelation_matrix_type;
           obj.user_speed                   = UE_params.user_speed;
           obj.realization_num              = UE_params.realization_num;
           obj.realization_num_total        = UE_params.realization_num_total;
           obj.CDD                          = UE_params.CDD;
           obj.mode                         = UE_params.mode;
           obj.nRX                          = UE_params.nRX;
           obj.timing_offset                = UE_params.timing_offset;
           obj.carrier_freq_offset          = UE_params.carrier_freq_offset;
           obj.timing_sync_method           = UE_params.timing_sync_method;
           obj.freq_sync_method             = UE_params.freq_sync_method;
           obj.rfo_correct_method           = UE_params.rfo_correct_method;
           obj.receiver                     = UE_params.receiver;
           obj.receiver_k                   = UE_params.receiver_k;
           if strcmp(fading,'BlockFading')
                obj.previous_channels            = zeros(Ntot,Nsub,10,UE_params.nRX,nTX);
           else
               obj.previous_channels            = zeros(Ntot,Nsub,1,UE_params.nRX,nTX);
           end
           obj.PMI_fb_gran                  = UE_params.PMI_fb_granularity;
           obj.CQI_fb_gran                  = UE_params.CQI_fb_granularity;
           obj.PMI_fb                       = UE_params.PMI_fb;
           obj.CQI_fb                       = UE_params.CQI_fb;
           obj.RIandPMI_fb                  = UE_params.RIandPMI_fb;
           obj.predict                      = UE_params.predict;
           obj.MSE                          = zeros(Ntot,Nsub);
           switch UE_params.SINR_averaging.averager
               case 'EESM'
                   obj.SINR_averager        = network_elements.eesmAverager(UE_params.SINR_averaging.EESMbetas,UE_params.SINR_averaging.MCSs);
               case 'MIESM'
                   obj.SINR_averager        = network_elements.miesmAverager(UE_params.SINR_averaging.MIESMbetas,UE_params.SINR_averaging.MCSs);
               otherwise
                   error('SINR averager not supported');
           end
           
           % Optional initialization 
           if isfield(UE_params,'clock')
               obj.clock                    = UE_params.clock;
           else
               obj.clock                    = 0;
           end
       end
  
       
   end
end 
