% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear 
clear global
%close all
clc

%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;

%% SNR setting
% SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
% SNR_stepsize = 1;
% SNR_window = 0.25;

Simulation_type = 'LTE_journal_paper';     %'SUSISO'
                                    %'MUSISO'
                                    %'SUMIMO'
                                    %'MUMIMO'
                                    %'SUSISO_quick_test'
                                    %'SUSISO_BLER_curves_batch'
                                    %'SUSISO_best_cqi'
                                    %'SUMIMO_quick_test'
                                    %'winner_model_example'
                                    %'wsa_2010_michal'
                                    
                                    

N_Ue = 2;
N_Bs = 2;
tx_mode = 3;
N_rx = 2;
N_tx = 2;
channel_type = 'PedA'
SNR_vec = 0:5:30;
Power_diff = [15;30;15;30];    % difference of the relativ signal power from one eNodeB to the other  [N_Ue*N_Bs x 1]
channel_est = 'PERFECT';
equalizer = 'SSD';
                    
%% Actual simulations
for cqi_i = 9

      N_subframes = 10;
      
        connection_table = false(N_Bs,N_Bs*N_Ue);
        connection_table(1,1:2) = 1;
        connection_table(2,3:4) = 1;
      

        LTE_load_parameters;  % Single User Multiple Input Multiple Output

        
    LTE_sim_main
    % Code to generate the output filename
%     output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
    output_filename = ['SISO_AWGN_CQI_feedback_max_TP'];
%     filename_suffix = [];

%     save(fullfile([output_filename '.mat']));
%     setpref('Internet','SMTP_Server','smtp.nt.tuwien.ac.at');
%     setpref('Internet','E_mail','stefan.schwarz@nt.tuwien.ac.at');
%     sendmail('stefan-schwarz@gmx.at','simulation','here it is',output_filename);
end
% shutdown(10)

