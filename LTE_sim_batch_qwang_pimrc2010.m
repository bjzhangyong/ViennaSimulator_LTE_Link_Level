% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear all
clear global
close all
clc

%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;

%% SNR setting

Simulation_type = 'pimrc_2010_qwang';     %'SUSISO'
                                    %'MUSISO'
                                    %'SUMIMO'
                                    %'MUMIMO'
                                    %'SUSISO_quick_test'
                                    %'SUSISO_BLER_curves_batch'
                                    %'SUSISO_best_cqi'
                                    %'SUMIMO_quick_test'
                                    %'winner_model_example'
                                    %'wsa_2010_michal'

%% simulation 1: over SNR
cqi_i = 9;
N_subframes = 1000;
channel_type_vec = {'AWGN','PedB'};
TX_mode_vec = {'SISO', 'MIMO'};
sync_mode_vec = {'perfect','compensated'};

line_style = {'b--','b-','g--','g-','r--','r-','m--','m-'};
sim_count = 1;

for channel_type_i = 1:length(channel_type_vec)
    channel_type = channel_type_vec{channel_type_i};
    
    for TX_mode_i = 1:length(TX_mode_vec)      
        TX_mode = TX_mode_vec{TX_mode_i};
        switch TX_mode
            case 'SISO'
                nTX = 1;
                nRX = 1;
                mode = 1;
            case 'MIMO'
                nTX = 4;
                nRX = 2;
                mode = 3;
            otherwise
                error('invalid TX mode.')
        end
        
        for sync_mode_i = 1:length(sync_mode_vec)
            sync_mode = sync_mode_vec{sync_mode_i};
            switch sync_mode
                case 'perfect'
                    CFO = 0;
                    freq_sync_method = 'perfect';
                case 'compensated'
                    CFO = pi;
                    freq_sync_method = 'estimated';
                otherwise
                    error('invalid frequency synchronization method.')
            end
            
            SNR_vec = 0:2:30;
            
            LTE_load_parameters;
            
            LTE_sim_main
            
            % result handling
%             output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
%             filename_suffix = [];
%             save(fullfile('./results',[output_filename filename_suffix '.mat']));

            if strcmp(sync_mode, 'compensated')
                ideal_int = round(CFO);
                ideal_frac = CFO - ideal_int;
                
                figure(1)
                semilogy(SNR_vec,mean((simulation_results.UE_specific.freq_offset_est_frac(1:5:end,:)-ideal_frac).^2,1),line_style{sim_count});
                hold on
                semilogy(SNR_vec,mean((simulation_results.UE_specific.freq_offset_est_frac+simulation_results.UE_specific.freq_offset_est_res-ideal_frac).^2,1),line_style{sim_count});
                xlabel('SNR [dB]')
                ylabel('Mean Squared Error')
                grid on
                
                figure(2)
                semilogy(SNR_vec,sum(logical(simulation_results.UE_specific.freq_offset_est_int(1:5:end,:)-ideal_int),1)/N_subframes*5,line_style{sim_count});
                hold on
                xlabel('SNR [dB]')
                ylabel('integer part error probability')
                grid on
            end
            
            figure(3)
            semilogy(SNR_vec,simulation_results.cell_specific.BER_uncoded_overall,line_style{sim_count})
            xlabel('SNR [dB]')
            ylabel('uncoded BER')
            hold on
            grid on
            
            figure(4)
            plot(SNR_vec,mean(sum(simulation_results.cell_specific.throughput_coded,3),1)/1e3,line_style{sim_count})
            xlabel('SNR [dB]')
            ylabel('throughput [Mbit/s]')
            hold on
            grid on
            
            sim_count = sim_count + 1;
        end
    end    
end


%% simulation 2: over CFO
N_subframes = 100;
channel_type_vec = {'AWGN','PedB'};
nTX = 1;
nRX = 1;
mode = 1;
freq_sync_method = 'none';
CFO_vec = logspace(-4,-1,100);

line_style = {'b-','r-'};
throughput_coded_cfo = zeros(length(CFO_vec),1);

for channel_type_i = 1:length(channel_type_vec)
    channel_type = channel_type_vec{channel_type_i};
    
    for cqi_i = 1:15
        
        for cfo_i = 1:length(CFO_vec)
            CFO = CFO_vec(cfo_i);
            
            SNR_vec = 30;
            
            LTE_load_parameters;
            
            LTE_sim_main
            
            throughput_coded_cfo(cfo_i) = mean(sum(simulation_results.cell_specific.throughput_coded,3),1)/1e3;
            
        end
        
        figure(sim_count)
        semilogx(CFO_vec,throughput_coded_cfo,line_style{channel_type_i})
        xlabel('SNR [dB]')
        ylabel('throughput [Mbit/s]')
        hold on
        grid on
        
    end
    sim_count = sim_count + 1;
end
% shutdown(10)

