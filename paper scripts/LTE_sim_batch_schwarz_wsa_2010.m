% Batch Simulation script to reproduce the results shown in 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mutual Information based Calculation of the
% Precoding Matrix Indicator for 3GPP UMTS/LTE 
% authors: Stefan Schwarz, Martin Wrulich, Markus Rupp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2,3,4,5

clear
clear global
close all
clc
cd ..
%% DEBUG level
global DEBUG_LEVEL;
DEBUG_LEVEL = 4; % 1-5

%% Actual simulations
Simulation_type = 'wsa_2010_schwarz';
receiver = 'ZF';
cqi_i = 4;
SNR_vec = [-4:12];
N_subframes = 5000;
for i1 = 1:2
    switch i1
        case 1
          for i2 = 1:3
              clearvars -except i1 N_subframes SNR_vec cqi_i DEBUG_LEVEL receiver i2 Simulation_type
              channel = 'PedA';
              nRX = 1;
              nTX = 2;
              i_max = 4;
              switch i2
                  case 1
                    PMI_fb = true;
                    LTE_load_parameters;
                    LTE_sim_main_single;
                    output_filename = 'PedA_2x1_feedback';
                    save(fullfile([output_filename '.mat']));
                  case 2
                    PMI_fb = false;
                    LTE_load_parameters;
                    LTE_sim_main_single;
                    output_filename = 'PedA_2x1_fixed';
                    save(fullfile([output_filename '.mat']));
                  case 3
                    PMI_fb = true;
                    LTE_load_parameters;
                    LTE_sim_main_single_optimumprecoder;
                    output_filename = 'PedA_2x1_optimal';
                    save(fullfile([output_filename '.mat']));
              end
          end
        case 2
          for i2 = 1:3
              clearvars -except i1 N_subframes SNR_vec cqi_i DEBUG_LEVEL receiver i2 Simulation_type
              channel = 'VehA';
              nRX = 2;
              nTX = 4;
              i_max = 16;
              switch i2
                  case 1
                    PMI_fb = true;
                    LTE_load_parameters;
                    LTE_sim_main_single;
                    output_filename = 'VehA_4x2_feedback';
                    save(fullfile([output_filename '.mat']));
                  case 2
                    PMI_fb = false;
                    LTE_load_parameters;
                    LTE_sim_main_single;
                    output_filename = 'VehA_4x2_fixed';
                    save(fullfile([output_filename '.mat']));
                  case 3
                    PMI_fb = true;
                    LTE_load_parameters;
                    LTE_sim_main_single_optimumprecoder;
                    output_filename = 'VehA_4x2_optimal';
                    save(fullfile([output_filename '.mat']));
              end
          end
    end
end

%% plot the results
% PedA
load('PedA_2x1_feedback.mat');
BERPedA_fb = simulation_results.cell_specific.BER_uncoded_overall;
BLERPedA_fb = simulation_results.cell_specific.BLER_overall;
load('PedA_2x1_fixed.mat');
BERPedA_fix = simulation_results.cell_specific.BER_uncoded_overall;
BLERPedA_fix = simulation_results.cell_specific.BLER_overall;
load('PedA_2x1_optimal.mat');
BERPedA_opt = simulation_results.cell_specific.BER_uncoded_overall;
BLERPedA_opt = simulation_results.cell_specific.BLER_overall;
figure(1)
clf;
semilogy(SNR_vec,BERPedA_fb,'r-x','linewidth',2,'Markersize',22);
hold on
grid on
semilogy(SNR_vec,BERPedA_fix,'k:x','linewidth',2,'Markersize',14);
semilogy(SNR_vec,BERPedA_opt,'k--+','linewidth',2,'Markersize',14);
set(gca,'Fontsize',24);
set(gca,'Linewidth',1.5);
xlab = xlabel('$\frac{E_s}{N_0}$ [dB]');
set(xlab,'interpreter','latex');
ylabel('BER');
legend('Mutual information feedback','Fixed feedback','Optimal choice');
title('2x1 PedA, ZF');

figure(2)
clf;
semilogy(SNR_vec,BLERPedA_fb,'r-x','linewidth',2,'Markersize',22);
hold on
grid on
semilogy(SNR_vec,BLERPedA_fix,'k:x','linewidth',2,'Markersize',14);
semilogy(SNR_vec,BLERPedA_opt,'k--+','linewidth',2,'Markersize',14);
set(gca,'Fontsize',24);
set(gca,'Linewidth',1.5);
xlab = xlabel('$\frac{E_s}{N_0}$ [dB]');
set(xlab,'interpreter','latex');
ylabel('BLER');
legend('Mutual information feedback','Fixed feedback','Optimal choice');
title('2x1 PedA, ZF');

%VehA
load('VehA_4x2_feedback.mat');
BERVehA_fb = simulation_results.cell_specific.BER_uncoded_overall;
BLERVehA_fb = simulation_results.cell_specific.BLER_overall;
load('VehA_4x2_fixed.mat');
BERVehA_fix = simulation_results.cell_specific.BER_uncoded_overall;
BLERVehA_fix = simulation_results.cell_specific.BLER_overall;
load('VehA_4x2_optimal.mat');
BERVehA_opt = simulation_results.cell_specific.BER_uncoded_overall;
BLERVehA_opt = simulation_results.cell_specific.BLER_overall;
figure(3)
clf;
semilogy(SNR_vec,BERVehA_fb,'r-x','linewidth',2,'Markersize',22);
hold on
grid on
semilogy(SNR_vec,BERVehA_fix,'k:x','linewidth',2,'Markersize',14);
semilogy(SNR_vec,BERVehA_opt,'k--+','linewidth',2,'Markersize',14);
set(gca,'Fontsize',24);
set(gca,'Linewidth',1.5);
xlab = xlabel('$\frac{E_s}{N_0}$ [dB]');
set(xlab,'interpreter','latex');
ylabel('BER');
legend('Mutual information feedback','Fixed feedback','Optimal choice');
title('4x2 VehA, ZF');

figure(4)
clf;
semilogy(SNR_vec,BLERVehA_fb,'r-x','linewidth',2,'Markersize',22);
hold on
grid on
semilogy(SNR_vec,BLERVehA_fix,'k:x','linewidth',2,'Markersize',14);
semilogy(SNR_vec,BLERVehA_opt,'k--+','linewidth',2,'Markersize',14);
set(gca,'Fontsize',24);
set(gca,'Linewidth',1.5);
xlab = xlabel('$\frac{E_s}{N_0}$ [dB]');
set(xlab,'interpreter','latex');
ylabel('BLER');
legend('Mutual information feedback','Fixed feedback','Optimal choice');
title('4x2 VehA, ZF');
% shutdown(10)

