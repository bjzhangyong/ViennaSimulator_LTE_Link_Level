% Reproducibility script for teh figure from the paper
% M. Simko, S. Pendl, S. Schwarz, Q. Wang, J. Colom Ikuno, M. Rupp: 
%"Optimal Pilot Symbol Power Allocation in LTE"; 
%Talk: 74th IEEE Vehicular Technology Conference (VTC2011-Fall), San Francisco, USA; 
%09-05-2011 - 09-08-2011; in: "Proc. 74th IEEE Vehicular Technology Conference (VTC2011-Fall)", (2011).

clear
clc
close all
%% Figure 2

p_off_dB = -10:0.1:10;
p_off = 10.^(-p_off_dB/20);

N_d = 960;
N_p = 48;
c_e = 0.3704;

f_o = ((p_off*N_d + N_p).^2) .* ((1./p_off.^2)+c_e);
figure
hold
plot(p_off_dB,f_o,'b')

N_d = 912;
N_p = 96;
c_e = 0.3704;

f_o = ((p_off*N_d + N_p).^2) .* ((1./p_off.^2)+c_e);
plot(p_off_dB,f_o,'r')

N_d = 864;
N_p = 144;
c_e = 0.5556;

f_o = ((p_off*N_d + N_p).^2) .* ((1./p_off.^2)+c_e);
plot(p_off_dB,f_o,'g')


N_d = 960;
N_p = 48;
c_e = 0.0394;

f_o = ((p_off*N_d + N_p).^2) .* ((1./p_off.^2)+c_e);
plot(p_off_dB,f_o,'Color','b','LineStyle','--')

N_d = 912;
N_p = 96;
c_e = 0.0394;

f_o = ((p_off*N_d + N_p).^2) .* ((1./p_off.^2)+c_e);
plot(p_off_dB,f_o,'Color','r','LineStyle','--')

N_d = 864;
N_p = 144;
c_e = 0.0544;

f_o = ((p_off*N_d + N_p).^2) .* ((1./p_off.^2)+c_e);
plot(p_off_dB,f_o,'Color','g','LineStyle','--')
axis([-10 10 1e6 1.8e6])
grid on
xlabel('p_off[dB]')
ylabel('f(o)');
%% Figure 3
cd ..
%path(cd,fileparts(pwd));
global DEBUG_LEVEL 
DEBUG_LEVEL = 4;
Simulation_type = 'LTE_journal_paper';
N_Ue = 1;
N_Bs = 1;
cqi_i = 9;
equalizer = 'ZF';
connection_table = true;
channel_type = 'PedA'
Power_diff = 100;
N_subframes = 500;
SNR_vec = [14];
p_off_dB = -10:0.5:10;
tho_matrix = nan(length(p_off_dB),N_subframes,3,2);% 3 correspodns to 3 different mimo modes, and 2 to 2 differetn channel estimators

for ce_i = 1:2
    
    if ce_i == 1
        % LS estimator
        channel_est = 'LS';
    else
        % LMMSE estimator
        channel_est = 'MMSE';
    end
    
    for mimo_i = 1:3
        if mimo_i == 1
            %SISO case
            tx_mode = 1;
            N_rx = 1;
            N_tx = 1;
        elseif mimo_i == 2
            %2x2 case
            tx_mode = 3;
            N_rx = 2;
            N_tx = 2;
        else
            %4x4 case
            tx_mode = 3;
            N_rx = 4;
            N_tx = 4;
        end
        
        for offset_i = 1:length(p_off_dB)
            pilot_offset = p_off_dB(offset_i);
        

            LTE_load_parameters;  % Single User Multiple Input Multiple Output
            LTE_sim_main
            tho_matrix(offset_i,:,mimo_i,ce_i) = sum(simulation_results.cell_specific.throughput_coded,3)/1000;
        end
        
    end
end

figure
plot(p_off_dB,mean(tho_matrix(:,:,1,1),2),'color','b','LineStyle','-')
hold on
plot(p_off_dB,mean(tho_matrix(:,:,2,1),2),'color','r','LineStyle','-')
plot(p_off_dB,mean(tho_matrix(:,:,3,1),2),'color','g','LineStyle','-')

plot(p_off_dB,mean(tho_matrix(:,:,1,2),2),'color','b','LineStyle','--')
plot(p_off_dB,mean(tho_matrix(:,:,2,2),2),'color','r','LineStyle','--')
plot(p_off_dB,mean(tho_matrix(:,:,3,2),2),'color','g','LineStyle','--')
grid on
xlabel('p_off[dB]')
ylabel('throughput [Mbit/s]');