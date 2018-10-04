function plot_quick_test_results_r553
% Script that shows example BLER and throughput results of the simulation.
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

folder = './examples/LTE_sim_batch_quick_test_r1089/';

flat_rayleigh_0retx = {
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode1_1x1_SSD_freqoff0_subframe_5000TTI_r1089_20110928_132404','SISO';
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode2_2x1_SSD_freqoff0_subframe_5000TTI_r1089_20110928_144708','TxD 2x1';
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode2_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20110928_161940','TxD 4x2';
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode3_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20110928_185415','OLSM 4x2';
    };

PedB_0retx = {
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode1_1x1_SSD_freqoff0_subframe_5000TTI_r1089_20110926_101006','SISO';
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode2_2x1_SSD_freqoff0_subframe_5000TTI_r1089_20110926_120642','TxD 2x1';
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode2_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20110926_141858','TxD 4x2';
    'cqi_7_0re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode3_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20110926_180648','OLSM 4x2';
    };

flat_rayleigh_3retx = {
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode1_1x1_SSD_freqoff0_subframe_5000TTI_r1089_20110928_192708','SISO';
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode2_2x1_SSD_freqoff0_subframe_5000TTI_r1089_20110928_204214','TxD 2x1';
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode2_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20110928_220506','TxD 4x2';
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_flat Rayleigh_TX_mode3_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20111004_135803','OLSM 4x2';
    };

PedB_3retx = {
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode1_1x1_SSD_freqoff0_subframe_5000TTI_r1089_20110926_185232','SISO';
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode2_2x1_SSD_freqoff0_subframe_5000TTI_r1089_20110926_203817','TxD 2x1';
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode2_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20110926_223655','TxD 4x2';
    'cqi_7_3re-tx_1TTI_UL_delay_1.4MHz_PedB_TX_mode3_4x2_SSD_freqoff0_subframe_5000TTI_r1089_20111004_143340','OLSM 4x2';
    };

plot_one_set(folder,flat_rayleigh_0retx,'flat rayleigh, 0 re-tx');
plot_one_set(folder,PedB_0retx,'PedB, 0 re-tx');
plot_one_set(folder,flat_rayleigh_3retx,'flat rayleigh, 3 re-tx');
plot_one_set(folder,PedB_3retx,'PedB, 3 re-tx');

function plot_one_set(folder,filenames,titles)

first_figure = figure;
first_axes = axes;
set(first_axes,'YScale','log');
title(first_axes,sprintf('BLER, CQI 7, PedB, 5000 subframes, %s',titles));
ylabel(first_axes,'BLER');
xlabel(first_axes,'SNR [dB]');
ylim([1e-3 1]);
hold all;
grid on;

second_figure = figure;
second_axes =  axes;
title(second_axes,sprintf('throughput, CQI 7, PedB, 5000 subframes, %s',titles));
ylabel(second_axes,'throughput [Mbps]');
xlabel(second_axes,'SNR [dB]');
hold all;
grid on;

% markers = {'x' '*' '+' '.'};
markers  = {'.' '.' '.' '.'};
marker_s = [ 15   15   15   15 ];

for file_idx=1:size(filenames,1)
    load(fullfile(folder,[ filenames{file_idx,1} '.mat']),'simulation_results','name','N_subframes');
    
    if length(markers)<file_idx
        current_marker = markers{end};
        current_marker_size = marker_s(end);
    else
        current_marker = markers{file_idx};
        current_marker_size = marker_s(file_idx);
    end
    displayname = filenames{file_idx,2};
    plot(first_axes,simulation_results.SNR_vector,simulation_results.cell_specific.BLER_overall,'DisplayName',displayname,'Marker',current_marker,'LineWidth',1,'MarkerSize',current_marker_size);
    plot(second_axes,simulation_results.SNR_vector,squeeze(sum(sum(simulation_results.cell_specific.throughput_coded,1),3))/(N_subframes*1e-3)/1e6,'DisplayName',displayname,'Marker',current_marker,'LineWidth',1,'MarkerSize',current_marker_size);
end

xlim(first_axes,[-10 20]);
xlim(second_axes,[-10 20]);
legend(first_axes,'show','Location','Best');
legend(second_axes,'show','Location','Best');
