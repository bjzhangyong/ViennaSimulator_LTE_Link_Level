% Batch file to reproduce the Figures shown in 
% "Throughput Maximizing Multiuser Scheduling with Adjustable Fairness" S. Schwarz, C. Mehlführer and M. Rupp

clear all
clear global
close all
clc
cd ..
current_dir = pwd;
Figures = 'Fig2'; % choose Fig2, Fig3 or Fig4
%% DEBUG level
global DEBUG_LEVEL LTE_params;
DEBUG_LEVEL = 4;

%% Actual simulations
cqi_i = 4;
N_subframes = 2500;  

%% Parameter definition
LTE_params.nBS = 1;  
LTE_params.introduce_frequency_offset =  false;
LTE_params.UE_config.channel_estimation_method = 'PERFECT';                       
LTE_params.UE_config.receiver = 'ZF';
LTE_params.UE_config.user_speed = 0/3.6;   
LTE_params.UE_config.carrier_freq_offset = pi;   
LTE_params.UE_config.freq_sync_method = 'perfect';
LTE_params.UE_config.rfo_correct_method = 'subframe'; 
LTE_params.ChanMod_config.filtering = 'BlockFading';  
LTE_params.scheduler.cqi  = 'set';
LTE_params.scheduler.PMI  = 2; 
LTE_params.uplink_delay = 0;
LTE_params.show_plots = false;
LTE_params.plot_confidence_intervals = true;
LTE_params.trace_subcarrier_SNR = false;
LTE_params.use_seed_reset = false; 
LTE_params.N_seed_reset = 1;     
LTE_params.carrier_freq = 2.1e9;   
LTE_params.speed_of_light = 299792458;
LTE_params.HARQ_processes = 8;           
LTE_params.max_HARQ_retransmissions = 0; 
LTE_params.SubcarrierSpacing = 15e3;
LTE_params.CyclicPrefix = 'normal';  
LTE_params.simulate_with_all_zero_sequences = false; 
LTE_params.random_noise_seeding = false;
LTE_params.noise_seed = 0;              
LTE_params.channel_matrix_source     = 'generated'; 
LTE_params.store_channel_trace       = false;         
LTE_params.channel_matrix_tracefile  = 'auto';       
LTE_params.CQI_mapping.coeffs = [0.5223 4.6176];    
LTE_params.CQI_mapping.table = [-500;-6.934;-5.147;-3.18;-1.254;0.761;2.70;4.697;6.528;8.576;10.37;12.3;14.18;15.89;17.82;19.83;21]; 
LTE_params.random_data_seeding = false;  
LTE_params.data_seed = 10;   
LTE_params.random_channel_param_seeding = false;  
LTE_params.channel_param_seed = 175;    
LTE_params.UE_config.LLR_clipping = 100;
LTE_params.UE_config.turbo_iterations = 8;                      
LTE_params.UE_config.N_soft = 1000000*LTE_params.HARQ_processes; 
LTE_params.UE_config.channel_interpolation_method = 'linear';    
LTE_params.UE_config.autocorrelation_matrix_type 	= 'ideal';   
LTE_params.UE_config.realization_num = 0;          
LTE_params.UE_config.realization_num_total = 20;   
LTE_params.UE_config.CDD =0;                       
LTE_params.UE_config.PMI_fb_granularity = 1;      
LTE_params.UE_config.CQI_fb_granularity = 1;       
LTE_params.UE_config.PMI_fb = true;             
LTE_params.UE_config.RIandPMI_fb = true;           
LTE_params.UE_config.CQI_fb = true;                 
LTE_params.UE_config.predict = false;               
LTE_params.UE_config.SINR_averaging.averager = 'MIESM';    
LTE_params.UE_config.SINR_averaging.EESMbetas = [5,5.01,5.01,0.84,1.67,1.61,1.64,3.87,5.06,6.4,12.59,17.59,23.33,29.45,33.05,35.41];  
LTE_params.UE_config.SINR_averaging.MIESMbetas = [4,3.07,4.41,0.6,1.16,1.06,1.06,0.87,1.01,1.04,1.03,1.11,1.01,1.07,1,1.05];  
LTE_params.UE_config.SINR_averaging.MCSs = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
LTE_params.UE_config.timing_offset = 23;   % timing offset in number of time samples
LTE_params.UE_config.timing_sync_method = 'perfect';% 'perfect','none', 'autocorrelation'
LTE_params.Pilots_power_offset = 0;
LTE_params.feedback.ignore_channel_estimation = false; 
LTE_params.feedback.channel_averaging = true; 
LTE_params.scheduler.av_window = 100;  
% LTE_params.scheduler.weights = 1/LTE_params.nUE*ones(LTE_params.nUE,1);
LTE_params.scheduler.stepsize = 0.5; 
LTE_params.ChanMod_config.interpolation_method = 'shift_to_nearest_neighbor'; 
LTE_params.ChanMod_config.corr_coefRX = 0.3;
LTE_params.ChanMod_config.corr_coefTX = 0.3;
LTE_params.ChanMod_config.sin_num = 10; 
LTE_params.eig_tresh_time = 0.1;
LTE_params.eig_tresh_freq = 0.1;
LTE_params.ChanMod_config.winner_settings.Scenario                 = 11;
LTE_params.ChanMod_config.winner_settings.PropagCondition          = 'NLOS';                           
LTE_params.ChanMod_config.winner_settings.SampleDensity            = 2;                                
LTE_params.ChanMod_config.winner_settings.UniformTimeSampling      = 'yes';                            
LTE_params.ChanMod_config.winner_settings.FixedPdpUsed             = 'no';            	               
LTE_params.ChanMod_config.winner_settings.FixedAnglesUsed          = 'no';                            
LTE_params.ChanMod_config.winner_settings.PolarisedArrays          = 'yes';                            
LTE_params.ChanMod_config.winner_settings.TimeEvolution            = 'no';                             
LTE_params.ChanMod_config.winner_settings.PathLossModelUsed        = 'no';                            
LTE_params.ChanMod_config.winner_settings.ShadowingModelUsed       = 'no';                             
LTE_params.ChanMod_config.winner_settings.PathLossModel            = 'pathloss';      	               
LTE_params.ChanMod_config.winner_settings.PathLossOption           = 'CR_light';                      
LTE_params.ChanMod_config.winner_settings.RandomSeed               = [];                              
LTE_params.ChanMod_config.winner_settings.UseManualPropCondition   = 'yes';                            
LTE_params.usePBCH = false;
LTE_params.usePDCCH = false;
LTE_params.trafficmodel.usetraffic_model = false; 
LTE_params.UE_config.mode = 1;
LTE_params.UE_config.nRX = 1;  
LTE_params.BS_config.nTx = 1;
LTE_params.simulation_type = 'normal';  
LTE_params.Bandwidth = 1.4e6;   
LTE_params.ChanMod_config.type = 'TU'; 
LTE_params.scheduler.PMI  = 2; 


LTE_params.Simulation_type = 'not_defined';
LTE_params.introduce_timing_offset = false;
LTE_params.SNR_estimation = false;

switch Figures
    case 'Fig2'
        SNR_vec = [15;0];
        LTE_params.nUE = 2; 
        LTE_params.ChanMod_config.time_correlation = 'independent';
        for loop1 = 1:7
            LTE_params.scheduler.PMI  = 2; 
            switch loop1
                case 1
                    LTE_params.scheduler.alpha = 7; 
                    LTE_params.scheduler.fairness = 0.95;
                    LTE_params.scheduler.type = 'best cqi';
                    LTE_params.scheduler.assignment = 'dynamic';
                    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'BCQI';
                    LTE_sim_main;
                    save(fullfile([output_filename '.mat']));
                case 2
                    LTE_params.scheduler.alpha = 7; 
                    LTE_params.scheduler.fairness = 0.95;
                    LTE_params.scheduler.type = 'max min';
                    LTE_params.scheduler.assignment = 'dynamic';
                    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'MM';
                    LTE_sim_main;
                    save(fullfile([output_filename '.mat']));
                case 3
                    LTE_params.scheduler.alpha = 7; 
                    LTE_params.scheduler.fairness = 0.95;
                    LTE_params.scheduler.type = 'proportional fair';
                    LTE_params.scheduler.assignment = 'dynamic';
                    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'PF';
                    LTE_sim_main;
                    save(fullfile([output_filename '.mat']));
                case 4
                    LTE_params.scheduler.alpha = 7; 
                    LTE_params.scheduler.fairness = 0.95;
                    LTE_params.scheduler.type = 'resource fair';
                    LTE_params.scheduler.assignment = 'dynamic';
                    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'RF';
                    LTE_sim_main;
                    save(fullfile([output_filename '.mat']));
                case 5
                    LTE_params.scheduler.alpha = 7; 
                    LTE_params.scheduler.fairness = 0.95;
                    LTE_params.scheduler.type = 'fixed';
                    LTE_params.scheduler.assignment = 'semi static';
                    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'RR';
                    LTE_sim_main;
                    save(fullfile([output_filename '.mat']));
                case 6
                    fairness = [0.5:0.025:0.975];
                    LTE_params.scheduler.type = 'var fair';
                    LTE_params.scheduler.assignment = 'dynamic';
                    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
                    LTE_params.scheduler.alpha = 7;
                    for loop2 = 1:length(fairness)
                        LTE_params.scheduler.PMI  = 2;
                        LTE_params.scheduler.fairness = fairness(loop2);
                        LTE_load_parameters_dependent;
                        LTE_load_parameters_generate_elements;
                        LTE_check_parameters;
                        output_filename = ['SOCP_', num2str(loop2)];
                        LTE_sim_main;
                        save(fullfile([output_filename '.mat']));
                    end
                case 7
                    alpha = [0,0.2,0.5,0.9,1.5,2,4,10,20,1000];
                    LTE_params.scheduler.type = 'var fair';
                    LTE_params.scheduler.assignment = 'dynamic';
                    LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
                    LTE_params.scheduler.fairness = 0.95;
                    for loop2 = 1:length(alpha)
                        LTE_params.scheduler.PMI  = 2;
                        LTE_params.scheduler.alpha = alpha(loop2);
                        LTE_load_parameters_dependent;
                        LTE_load_parameters_generate_elements;
                        LTE_check_parameters;
                        output_filename = ['Alpha_', num2str(loop2)];
                        LTE_sim_main;
                        save(fullfile([output_filename '.mat']));
                    end
            end
        end
        cd(current_dir);
        Schwarz_ICC_plot_fair;
    case 'Fig3'
        SNR_vec = [15;12;10;5;0];
        LTE_params.nUE = 5;  
        LTE_params.ChanMod_config.time_correlation = 'correlated'; 
        LTE_params.scheduler.type = 'variable fair';
        LTE_params.scheduler.assignment = 'dynamic';
        LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
        LTE_params.scheduler.fairness = 0.95; 
        LTE_params.scheduler.alpha = 7; 
        LTE_load_parameters_dependent;
        LTE_load_parameters_generate_elements;
        LTE_check_parameters;
        output_filename = 'Temp_correlated';
        LTE_sim_main;
        save(fullfile([output_filename '.mat']));
        cd(current_dir);
        Schwarz_ICC_plot_temp;
    case 'Fig4'
        SNR_vec = [15;12;10;5;0];
        LTE_params.nUE = 5;  
        LTE_params.ChanMod_config.time_correlation = 'independent'; 
        LTE_params.scheduler.type = 'variable fair';
        LTE_params.scheduler.assignment = 'dynamic';
        LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);
        LTE_params.scheduler.fairness = 0.95; 
        LTE_params.scheduler.alpha = 7; 
        LTE_load_parameters_dependent;
        LTE_load_parameters_generate_elements;
        LTE_check_parameters;
        output_filename = 'Temp_independent';
        LTE_sim_main;
        save(fullfile([output_filename '.mat']));
        cd(current_dir);
        Schwarz_ICC_plot_temp;
end







