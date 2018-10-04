 % Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear
%clear global
%close all
clc

%% DEBUG level
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;

%% SNR setting
% SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
% SNR_stepsize = 1;
% SNR_window = 0.25;
% power = [];
% noise = [];
Simulation_type = 'MUMIMO';     %'SUSISO'
                                    %'MUSISO'
                                    %'SUMIMO'
                                    %'MUMIMO'
                                    %'SUSISO_quick_test'
                                    %'SUSISO_BLER_curves_batch'
                                    %'SUSISO_best_cqi'
                                    %'SUMIMO_quick_test'
                                    %'winner_model_example'
                                    %'wsa_2010_michal'
                                    
                                    
%  counti = 1;
%  channel_estimation_error_freq_depend = zeros(72,14,1,2,500);
% Hsave = zeros(72,1000,32);
%% Actual simulations
for cqi_i = 9
%     N_subframes = [200;ones(200,1);200];    %     SNR_vec = SNR_30percent(cqi_i)-SNR_window*2.5:SNR_stepsize:SNR_30percent(cqi_i)+SNR_window;SNR_vec = 50;
      N_subframes = 100;
%     channel_estimation_error_freq_depend = zeros(72,14,1,2,N_subframes);   
%     SNR_vec = [7.7;5.72;27.63;14.65;12.37;8.12;15.66;26.27;8.05;28.31;15.41;15.84;15.85;8.81;-0.42;26.45;12.09;10.88;32.96;11.22];
SNR_vec = [10];
%       SNR_vec(:,:,1) = [0:10:30];
%       SNR_vec(:,:,1) = [15;0];
%       SNR_vec(:,:,3) = [15;0];
%       SNR_vec(:,:,2:201) = [linspace(15,5,200);linspace(0,0,200)];
%       SNR_vec(:,:,2:201) = [linspace(15,5,200);linspace(12,12,200);linspace(10,10,200);linspace(5,5,200);linspace(0,0,200)];
%       SNR_vec(:,:,202) = [5;12;10;5;0];
%       SNR_vec(:,:,202) = [5;0];
%       SNR_vec(:,:,103) = [15;0];
      LTE_load_parameters;  % Single User Multiple Input Multiple Output
%     LTE_params_36101;  % Multi User Multiple Input Multiple Output
%     LTE_load_parameters_SUSISO;  % Single User Single Input Single Output
%     LTE_load_parameters_MUMIMO_old;   % Multi User Single Input Single
% %     Output
%     i_max = 3;
%     LTE_sim_main_single_optimumprecoder
%     LTE_sim_main_par_36101;
%     LTE_sim_main_par_opt_RI_PMI_CQI
    LTE_sim_main
    % Code to generate the output filename
%     output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
    output_filename = ['Mobilkom_CLSM_2x2_xx_1/2'];
%     filename_suffix = [];

%     save(fullfile([output_filename '.mat']));
%     setpref('Internet','SMTP_Server','smtp.nt.tuwien.ac.at');
%     setpref('Internet','E_mail','stefan.schwarz@nt.tuwien.ac.at');
%     sendmail('stefan-schwarz@gmx.at','simulation','here it is',output_filename);
end
% shutdown(10)

