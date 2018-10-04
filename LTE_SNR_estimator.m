function [S_plus_noise_power N_power] = LTE_SNR_estimator(LTE_params,rx_signal);
%
%
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) by INTHFT
% www.nt.tuwien.ac.at
%
% input :
% output:   tx_streams                      ... struct containing TX data
% streams
%
% date of creation: 2011/04/04
% last changes: 2011/04/04  Simko



% noise power estimation
noise_start = abs(rx_signal(1:LTE_params.number_of_zeros_for_SNR_estimation,:)).^2;
noise_end = abs(rx_signal(end-LTE_params.number_of_zeros_for_SNR_estimation+1:end,:)).^2;
%nosise_power_vec = [noise_start; noise_end];
N_power = (1/(LTE_params.Nfft/LTE_params.Ntot))*[mean(noise_start,1); mean(noise_end,1)];


% signal power estimation
start_pos = ceil(length(rx_signal)/2);
nr_of_samples = floor(start_pos/LTE_params.number_of_zeros_for_SNR_estimation)-2;
S_plus_noise_power = mean(abs(rx_signal(start_pos-nr_of_samples*LTE_params.number_of_zeros_for_SNR_estimation:start_pos+nr_of_samples*LTE_params.number_of_zeros_for_SNR_estimation,:)).^2);