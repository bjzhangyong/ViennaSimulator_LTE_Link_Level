% Batch file to reproduce the Figures shown in 
% "Low Complexity Approximate Maximum Throughput Scheduling for LTE" S. Schwarz, C. Mehlführer and M. Rupp

clear all
clear global
close all
clc
cd ..
%% DEBUG level
global DEBUG_LEVEL LTE_params;
DEBUG_LEVEL = 4;

Figure_choice = 'Fig2_3'; % change this to the Figures from the paper you want to reproduce (Fig2_3 or Fig4_6)

%% Actual simulations
cqi_i = 4;
N_subframes = 1000;  

%% Parameter definition
LTE_params.nBS = 1;  
LTE_params.introduce_frequency_offset =  false;
LTE_params.UE_config.channel_estimation_method = 'PERFECT';     
LTE_params.UE_config.mode = 1;                    
LTE_params.UE_config.nRX = 1;                      
LTE_params.UE_config.receiver = 'ZF';
LTE_params.UE_config.user_speed = 0/3.6;   
LTE_params.UE_config.carrier_freq_offset = pi;   
LTE_params.UE_config.freq_sync_method = 'perfect';
LTE_params.UE_config.rfo_correct_method = 'subframe'; 
LTE_params.BS_config.nTx = 1;
LTE_params.ChanMod_config.filtering = 'BlockFading';  
LTE_params.scheduler.assignment = 'dynamic';
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
LTE_params.scheduler.fairness = 0.95; 
LTE_params.scheduler.alpha = 7;   
% LTE_params.scheduler.weights = 1/LTE_params.nUE*ones(LTE_params.nUE,1);
LTE_params.scheduler.stepsize = 0.5; 
LTE_params.ChanMod_config.interpolation_method = 'shift_to_nearest_neighbor'; 
LTE_params.ChanMod_config.corr_coefRX = 0.3;
LTE_params.ChanMod_config.corr_coefTX = 0.3;
LTE_params.ChanMod_config.sin_num = 10; 
LTE_params.ChanMod_config.time_correlation = 'independent';   
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

LTE_params.simulation_type = 'parallel'; 

LTE_params.Simulation_type = 'not_defined';
LTE_params.introduce_timing_offset = false;
LTE_params.SNR_estimation = false;
%%
switch Figure_choice
    case 'Fig2_3'
        LTE_params.nUE = 2;     
        LTE_params.Bandwidth = 1.4e6;   
        LTE_params.ChanMod_config.type = 'VehA'; 
        loop_nr = 3;
    case 'Fig4_6'
        LTE_params.nUE = 25;     
        LTE_params.Bandwidth = 10e6;   
        LTE_params.ChanMod_config.type = 'TU'; 
        loop_nr = 7;
    otherwise
        error('Figure choice not supported');
end
LTE_params.connection_table = true(LTE_params.nBS,LTE_params.nBS*LTE_params.nUE);

for loop_i = 1:loop_nr
    LTE_params.scheduler.PMI  = 2; 
    switch Figure_choice
        case 'Fig2_3'
            SNR_vec = [-10:2:30;-13:2:27];
            switch loop_i
                case 1
                    LTE_params.scheduler.type = 'best cqi';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig2_3_best_cqi';
                case 2
                    LTE_params.scheduler.type = 'max throughput';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig2_3_AMT';
                case 3
                    LTE_params.scheduler.type = 'Kwan';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig2_3_KMT';
            end
        case 'Fig4_6'
            SNR_vec = [1:25]';
            switch loop_i
                case 1
                    LTE_params.scheduler.type = 'best cqi';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig4_6_best_cqi';
                case 2
                    LTE_params.scheduler.type = 'max throughput';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig4_6_AMT';
                case 3
                    LTE_params.scheduler.type = 'Kwan';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig4_6_KMT';
                case 4
                    LTE_params.scheduler.type = 'proportional fair';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig4_6_PF';
                case 5
                    LTE_params.scheduler.type = 'resource fair';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig4_6_RF';
                case 6
                    LTE_params.scheduler.type = 'max min';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig4_6_MM';
                case 7
                    LTE_params.scheduler.type = 'fixed';
                    LTE_params.scheduler.assignment = 'semi static';
                    LTE_load_parameters_dependent;
                    LTE_load_parameters_generate_elements;
                    LTE_check_parameters;
                    output_filename = 'Fig4_6_RR';
            end
    end
    LTE_sim_main;
    save(fullfile([output_filename '.mat']));
end

switch Figure_choice
    case 'Fig2_3'
        fig = figure(1);
        set(fig,'units','normalized','outerposition',[0 0 1 1]);
        load('Fig2_3_best_cqi.mat')
        y_dat(1,:) = mean(sum(simulation_results.UE_specific(1).throughput_coded,3)/LTE_params.Tsubframe/1e6);
        y_dat(2,:) = mean(sum(simulation_results.UE_specific(2).throughput_coded,3)/LTE_params.Tsubframe/1e6);
        y_dat(3,:) = y_dat(2,:)+y_dat(1,:);
        plot_multiplex(SNR_vec,y_dat(3,:),'SNR [dB]','Throughput [Mbit/s]','Best CQI scheduler throughput, 1x1 VehA, 1.4MHz',1)
        hold on
        grid on
        
        fig = figure(2);
        set(fig,'units','normalized','outerposition',[0 0 1 1]);
        y_dat(1,:) = simulation_results.UE_specific(1).BLER_overall;
        y_dat(2,:) = simulation_results.UE_specific(2).BLER_overall;
        semilogy(SNR_vec(1,:),y_dat(1,:),'rx-','Linewidth',2.5,'Markersize',20)
        hold on
        grid on
        set(gca,'Fontsize',30)
        set(gca,'Linewidth',1.5)
        semilogy(SNR_vec(2,:),y_dat(2,:),'rx--','Linewidth',2.5,'Markersize',20)

        figure(1);
        clear y_dat
        load('Fig2_3_AMT.mat')
        y_dat(1,:) = mean(sum(simulation_results.UE_specific(1).throughput_coded(:,:,:),3)/LTE_params.Tsubframe/1e6);
        y_dat(2,:) = mean(sum(simulation_results.UE_specific(2).throughput_coded(:,:,:),3)/LTE_params.Tsubframe/1e6);
        y_dat(3,:) = y_dat(2,:)+y_dat(1,:);
        plot(SNR_vec(1,:),y_dat(3,:),'bd:','Linewidth',2.5,'Markersize',16);
        
        figure(2)
        clear y_dat
        y_dat(1,:) = simulation_results.UE_specific(1).BLER_overall;
        y_dat(2,:) = simulation_results.UE_specific(2).BLER_overall;
        semilogy(SNR_vec(1,:),y_dat(1,:),'bd-','Linewidth',2.5,'Markersize',16)
        semilogy(SNR_vec(2,:),y_dat(2,:),'bd--','Linewidth',2.5,'Markersize',16)
        
        figure(1);
        clear y_dat
        load('Fig2_3_KMT.mat')
        y_dat(1,:) = mean(sum(simulation_results.UE_specific(1).throughput_coded(:,:,:),3)/LTE_params.Tsubframe/1e6);
        y_dat(2,:) = mean(sum(simulation_results.UE_specific(2).throughput_coded(:,:,:),3)/LTE_params.Tsubframe/1e6);
        y_dat(3,:) = y_dat(2,:)+y_dat(1,:);
        plot(SNR_vec(1,:),y_dat(3,:),'go-','Linewidth',2.5,'Markersize',16);
        legend('BCQI','AMT','KMT','Location','Northwest');
        ylim([0,5.5])
        
        figure(2)
        y_dat(1,:) = simulation_results.UE_specific(1).BLER_overall;
        y_dat(2,:) = simulation_results.UE_specific(2).BLER_overall;
        semilogy(SNR_vec(1,:),y_dat(1,:),'go-','Linewidth',2.5,'Markersize',16)
        semilogy(SNR_vec(2,:),y_dat(2,:),'go--','Linewidth',2.5,'Markersize',16)
        xlim([-15,30]);
        xlabel('SNR [dB]');
        ylabel('BLER');
        legend('UE1 BCQI','UE2 BCQI','UE1 AMT','UE2 AMT','UE1 KMT','UE2 KMT','Location','Southwest');
        
    case 'Fig4_6'
        count = 0;

        load('Fig4_6_RR.mat');
        count = count+1;
        for ue_i = 1:25+1
            if ue_i < 26
                Y(ue_i,count) = mean(sum(simulation_results.UE_specific(ue_i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
            else
                Y(ue_i,count) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            end
        end
        load('Fig4_6_MM.mat');
        count = count+1;
        for ue_i = 1:25+1
            if ue_i < 26
                Y(ue_i,count) = mean(sum(simulation_results.UE_specific(ue_i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
            else
                Y(ue_i,count) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            end
        end
        load('Fig4_6_PF.mat');
        count = count+1;
        for ue_i = 1:25+1
            if ue_i < 26
                Y(ue_i,count) = mean(sum(simulation_results.UE_specific(ue_i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
            else
                Y(ue_i,count) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            end
        end
        load('Fig4_6_RF.mat');
        count = count+1;
        for ue_i = 1:25+1
            if ue_i < 26
                Y(ue_i,count) = mean(sum(simulation_results.UE_specific(ue_i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
            else
                Y(ue_i,count) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            end
        end
        load('Fig4_6_AMT.mat');
        count = count+1;
        for ue_i = 1:25+1
            if ue_i < 26
                Y(ue_i,count) = mean(sum(simulation_results.UE_specific(ue_i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
            else
                Y(ue_i,count) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            end
        end
        load('Fig4_6_best_cqi.mat');
        count = count+1;
        for ue_i = 1:25+1
            if ue_i < 26
                Y(ue_i,count) = mean(sum(simulation_results.UE_specific(ue_i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
            else
                Y(ue_i,count) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            end
        end
        load('Fig4_6_KMT.mat');
        count = count+1;
        for ue_i = 1:25+1
            if ue_i < 26
                Y(ue_i,count) = mean(sum(simulation_results.UE_specific(ue_i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
            else
                Y(ue_i,count) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            end
        end
        bar3(SNR_vec.',Y(1:end-1,:),0.5,'Detached');
        grid on
        xlim([1,7]);
        ylim([0,26]);
        zlim([0,10]);
        set(gca,'Fontsize',30);
        set(gca,'Linewidth',1.5);
        legend('RR','MaxMin','PF','RF','AMT','BCQI','KMT','Location','Northwest')
        ax1 = get(gcf);
        colors = ax1.Colormap;
        zlabel('Throughput [Mbit/s]');
        ylabel('UE index (= UE SNR [dB])');
        set(gca,'XTickLabel',[]);

        figure(2)
        set(gca,'Fontsize',30);
        set(gca,'Linewidth',1.5);
        for i = 1:size(Y,2)
            bar(i,Y(end,i),'FaceColor',colors(min(floor((i-1)*length(colors)/(size(Y,2)-1)+1),length(colors)),:))
            hold on
        end
        grid on
        set(gca,'XTick',1:size(Y,2))
        set(gca,'XTickLabel',{'RR','Max.Min.','PF','RF','AMT','BCQI','KMT'});
        xlabel('Scheduler');
        ylabel('Throughput [Mbit/s]');
        xlim([0,8])
        J = sum(Y(1:end-1,:),1).^2./(size(Y(1:end-1,:),1)*sum(Y(1:end-1,:).^2,1));
        figure(3)
        set(gca,'Fontsize',30);
        set(gca,'Linewidth',1.5);
        for i = 1:size(Y,2)
            bar(i,J(i),'FaceColor',colors(min(floor((i-1)*length(colors)/(size(Y,2)-1)+1),length(colors)),:))
            hold on
        end
        grid on
        xlim([0,8])
        set(gca,'XTick',1:size(Y,2))
        set(gca,'XTickLabel',{'RR','Max.Min.','PF','RF','AMT','BCQI','KMT'});
        xlabel('Scheduler');
        ylabel('Jain''s Fairness Index');
end

