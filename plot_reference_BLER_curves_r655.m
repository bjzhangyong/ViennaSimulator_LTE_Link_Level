% Script that shows reference BLER curves.
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

first_figure = figure;
first_axes = axes;
set(first_axes,'YScale','log');
title('BLER, 1.4MHz, SISO AWGN, 5000 subframes');
ylabel('BLER');
xlabel('SNR [dB]');
ylim([1e-3 1]);
hold all;
grid on;

second_figure = figure;
second_axes =  axes;
title(second_axes,sprintf('throughput, 1.4MHz, SISO AWGN, 5000 subframes'));
ylabel(second_axes,'throughput [Mbps]');
xlabel(second_axes,'SNR [dB]');
hold all;
grid on;

files = {
    './examples/AWGN_1.4MHz_r655/cqi_1_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100614_152540.mat'
    './examples/AWGN_1.4MHz_r655/cqi_2_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100614_192832.mat'
    './examples/AWGN_1.4MHz_r655/cqi_3_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100614_233339.mat'
    './examples/AWGN_1.4MHz_r655/cqi_4_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100615_034238.mat'
    './examples/AWGN_1.4MHz_r655/cqi_5_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100615_075750.mat'
    './examples/AWGN_1.4MHz_r655/cqi_6_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100615_155128.mat'
    './examples/AWGN_1.4MHz_r655/cqi_7_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100615_203112.mat'
    './examples/AWGN_1.4MHz_r655/cqi_8_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100616_011251.mat'
    './examples/AWGN_1.4MHz_r655/cqi_9_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100616_060628.mat'
    './examples/AWGN_1.4MHz_r655/cqi_10_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100616_112217.mat'
    './examples/AWGN_1.4MHz_r655/cqi_11_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100616_164805.mat'
    './examples/AWGN_1.4MHz_r655/cqi_12_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100616_222359.mat'
    './examples/AWGN_1.4MHz_r655/cqi_13_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100617_041145.mat'
    './examples/AWGN_1.4MHz_r655/cqi_14_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100617_100858.mat'
    './examples/AWGN_1.4MHz_r655/cqi_15_0re-tx_0TTI_UL_delay_1.4MHz_AWGN_TX_mode1_1x1_SSD_freqoff0_none_5000TTI_r655_20100617_161416.mat'
    };

for cqi_value=1:15
    load(files{cqi_value},'simulation_results','cqi_i','N_subframes');

    name = sprintf('CQI %02d',cqi_i);
    
    plot(first_axes,simulation_results.SNR_vector,simulation_results.cell_specific.BLER_overall,'DisplayName',name);
    plot(second_axes,simulation_results.SNR_vector,squeeze(sum(sum(simulation_results.cell_specific.throughput_coded,1),3))/(N_subframes*1e-3)/1e6,'DisplayName',name);
end

%xlim(first_axes,[-15 25]);
%xlim(second_axes,[-15 25]);
xlim(first_axes,'auto');
xlim(second_axes,'auto');
legend(first_axes,'show','Location','SouthEastOutside');
legend(second_axes,'show','Location','SouthEastOutside');
