% Basic batch simulation script
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

clear
clear global
% close all
clc
cd ..
%% DEBUG level
global DEBUG_LEVEL;
DEBUG_LEVEL = 4; % 1-5

%% SNR setting
SNR_30percent = [-7, -5, -3, -1, 1, 3, 3, 7, 9, 11, 13, 14.5, 16, 17.75, 19.5];
SNR_stepsize = 1;
SNR_window = 0.25;



N_subframes = 1000;
cqi_i = 10;

Simulation_type = 'wsa_2010_michal';     %'SUSISO'
                                    %'MUSISO'
                                    %'SUMIMO'
                                    %'MUMIMO'
                                    %'SUSISO_quick_test'
                                    %'SUSISO_BLER_curves_batch'
                                    %'SUSISO_best_cqi'
                                    %'SUMIMO_quick_test'
                                    %'winner_model_example'
                                    %'wsa_2010_michal'
                                    
%% Simulation over velocity with fixed SNR
speed_vec = [0:25:350]/3.6;
estimator_vec = 1:5;
tho = nan(length(speed_vec),length(estimator_vec));
mse = nan(length(speed_vec),length(estimator_vec));

for est_i = 1:length(estimator_vec)
    for speed_i = 1:length(speed_vec)
        speed = speed_vec(speed_i)
        
            
            SNR_vec = 20;
            switch estimator_vec(est_i)
                case 1
                    estimator = 'PERFECT'
                case 2
                    estimator = 'MMSE'
                case 3
                    estimator = 'ALMMSE'
                case 4
                    estimator = 'LS'
                case 5
                    estimator = 'LS_block'
            end
            LTE_load_parameters;  % Multi User Multiple Input Multiple Output using Winner II + Channel Model
            LTE_params.eig_tresh_time = 0.1;
            LTE_params.eig_tresh_freq = 0.1;

            % LTE_load_parameters_MUMIMO_winner_example;

            LTE_sim_main;
            tho(speed_i,est_i) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            mse(speed_i,est_i) = simulation_results.cell_specific.MSE_overall;

            % Code to generate the output filename
            output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
            filename_suffix = [];

        %     save(fullfile('./results',[output_filename filename_suffix '.mat']));
        
    end
end
figure
plot(speed_vec*3.6,tho(:,1),'Color','black','Marker','d');
hold on
plot(speed_vec*3.6,tho(:,2),'Color','b','Marker','+');
plot(speed_vec*3.6,tho(:,3),'Color','r','Marker','x');
plot(speed_vec*3.6,tho(:,4),'Color','g','Marker','s');
plot(speed_vec*3.6,tho(:,5),'Color','m','Marker','o');
xlabel('User speed [km/h]');
ylabel('Throughput [Mbit/s]');
legend('PERFECT','LMMSE','ALMMSE','LS','LS block')
axis([0 350 0 9])


figure
semilogy(speed_vec*3.6,mse(:,2),'Color','b','Marker','+');
hold on

semilogy(speed_vec*3.6,mse(:,3),'Color','r','Marker','x');
semilogy(speed_vec*3.6,mse(:,4),'Color','g','Marker','s');
semilogy(speed_vec*3.6,mse(:,5),'Color','m','Marker','o');
xlabel('User speed [km/h]');
ylabel('Throughput [Mbit/s]');
legend('LMMSE','ALMMSE','LS','LS block')

save(fullfile('./results',[output_filename filename_suffix '.mat']));

%% Simulation over SNR with fixed velocity
SNR_vec = 0:30;
tho_mat = nan(length(SNR_vec),length(estimator_vec));
mse_mat = nan(length(SNR_vec),length(estimator_vec));
for est_i = 1:length(estimator_vec)
    
    speed = 60/3.6;
    
        
        
        switch estimator_vec(est_i)
            case 1
                estimator = 'PERFECT'
            case 2
                estimator = 'MMSE'
            case 3
                estimator = 'ALMMSE'
            case 4
                estimator = 'LS'
            case 5
                estimator = 'LS_block'
        end
        LTE_load_parameters_MUMIMO_michal_wsa_2010;  % Multi User Multiple Input Multiple Output using Winner II + Channel Model
        LTE_params.simulation_type = 'parallel'; % 'parallel'
        LTE_params.eig_tresh_time = 0.1;
        LTE_params.eig_tresh_freq = 0.1;
        
        % LTE_load_parameters_MUMIMO_winner_example;
        
        LTE_sim_main;
        tho_mat(:,est_i) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
        mse_mat(:,est_i) = simulation_results.cell_specific.MSE_overall;
        
        % Code to generate the output filename
        output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
        filename_suffix = [];
        
        save(fullfile('./results',[output_filename filename_suffix '.mat']));
    
    
end

figure
plot(SNR_vec,tho_mat(:,1),'Color','black','Marker','d');
hold on
plot(SNR_vec,tho_mat(:,2),'Color','b','Marker','+');
plot(SNR_vec,tho_mat(:,3),'Color','r','Marker','x');
plot(SNR_vec,tho_mat(:,4),'Color','g','Marker','s');
plot(SNR_vec,tho_mat(:,5),'Color','m','Marker','o');
xlabel('SNR [dB]');
ylabel('Throughput [Mbit/s]');
legend('PERFECT','LMMSE','ALMMSE','LS','LS block')
axis([0 30 0 9])


figure
semilogy(SNR_vec,mse_mat(:,2),'Color','b','Marker','+');
hold on

semilogy(SNR_vec,mse_mat(:,3),'Color','r','Marker','x');
semilogy(SNR_vec,mse_mat(:,4),'Color','g','Marker','s');
semilogy(SNR_vec,mse_mat(:,5),'Color','m','Marker','o');
xlabel('SNR [dB]');
ylabel('MSE');
legend('LMMSE','ALMMSE','LS','LS block')

output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
filename_suffix = [];

save(fullfile('./results',[output_filename filename_suffix '.mat']));



%% Simulation of the ALMMSE with different number of eingevalues
speed_vec = [0:25:200]/3.6;
estimator = 'ALMMSE';  %estimator set to ALMMSE

%eigenvalue_tresh = [0.05,0.05,0.05,0.05,0.05,0.1,0.1,0.1,0.1,0.1,3,3,3,3,3;0.1,0.3,2,16,17,0.1,0.3,2,16,17,0.1,0.3,2,16,17,];
%eigenvalue_tresh = [0.05,0.05,0.05,0.1,3;
%                    0.1,0.3,2,0.1,17];
                    
eigenvalue_tresh = [0.05,0.1,0.1,0.1,0.1,3;
                    0.1,0.1,0.3,2,17,17];

tho_mat2 = nan(length(speed_vec),size(eigenvalue_tresh,2));
mse_mat2 = nan(length(speed_vec),size(eigenvalue_tresh,2));

for eig_i = 1:size(eigenvalue_tresh,2)
    eig_i
    for speed_i = 1:length(speed_vec)
        speed = speed_vec(speed_i)
        
            
            SNR_vec = 20;
            LTE_load_parameters_MUMIMO_michal_wsa_2010;  % Multi User Multiple Input Multiple Output using Winner II + Channel Model
            LTE_params.eig_tresh_time = eigenvalue_tresh(1,eig_i);
            LTE_params.eig_tresh_freq = eigenvalue_tresh(2,eig_i);
%             if speed == 25/3.6 && (eig_i == 2 || eig_i == 3 || eig_i == 4 || eig_i == 5 )
%                 LTE_params.eig_tresh_time = 0.05;
%             end

            % LTE_load_parameters_MUMIMO_winner_example;

            LTE_sim_main;
            tho_mat2(speed_i,eig_i) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
            mse_mat2(speed_i,eig_i) = simulation_results.cell_specific.MSE_overall;

            % Code to generate the output filename
            output_filename = LTE_common_generate_output_filename(LTE_params,N_subframes);
            filename_suffix = [];

        %     save(fullfile('./results',[output_filename filename_suffix '.mat']));
        
    end
end

figure
plot(speed_vec*3.6,tho_mat2(:,2),'Color','b','Marker','+');
hold on
plot(speed_vec*3.6,tho_mat2(:,3),'Color','r','Marker','x');
plot(speed_vec*3.6,tho_mat2(:,4),'Color','g','Marker','s');
plot(speed_vec*3.6,tho_mat2(:,6),'Color','y','Marker','h');
legend('N_freq = 5 N_time = 2','N_freq = 4 N_time = 2','N_freq = 3 N_time = 2','N_freq = 1 N_time = 1');
xlabel('User speed [km/h]');
ylabel('Throughput [Mbit/s]');
%LMMSE and LS
plot(speed_vec*3.6,tho(1:length(speed_vec),2),'Color','black','LineStyle',':');
plot(speed_vec*3.6,tho(1:length(speed_vec),4),'Color','black','LineStyle',':');

figure
semilogy(speed_vec*3.6,mse_mat2(:,2),'Color','b','Marker','+');
hold on
semilogy(speed_vec*3.6,mse_mat2(:,3),'Color','r','Marker','x');
semilogy(speed_vec*3.6,mse_mat2(:,4),'Color','g','Marker','s');
semilogy(speed_vec*3.6,mse_mat2(:,6),'Color','y','Marker','h');
legend('N_freq = 5 N_time = 2','N_freq = 4 N_time = 2','N_freq = 3 N_time = 2','N_freq = 1 N_time = 1');
xlabel('User speed [km/h]');
ylabel('MSE');
%LMMSE and LS
semilogy(speed_vec*3.6,mse(1:length(speed_vec),2),'Color','black','LineStyle',':');
semilogy(speed_vec*3.6,mse(1:length(speed_vec),4),'Color','black','LineStyle',':');

% figure
% color_mat = colormap(jet(size(eigenvalue_tresh,2)));
% hold on
% for fig_i = 1:size(eigenvalue_tresh,2)
%     plot(speed_vec*3.6,tho_mat2(:,fig_i),'color',color_mat(fig_i,:))
% end



% figure
% hold on
% for fig_i = 1:size(eigenvalue_tresh,2)
%     semilogy(speed_vec*3.6,mse_mat2(:,fig_i),'color',color_mat(fig_i,:))
% end


save(fullfile('./results',[output_filename filename_suffix '.mat']));