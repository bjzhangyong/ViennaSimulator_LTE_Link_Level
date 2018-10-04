function results = LTE_set_results(UE_ACK,UE_rv_idx,UE_RBs_assigned,UE_biterrors_coded,UE_biterrors_uncoded,UE_blocksize_coded,UE_blocksize_uncoded,...
    UE_throughput_coded,UE_throughput_uncoded,UE_throughput_useful,UE_FER_coded,UE_FER_uncoded,UE_used_codewords,UE_used_CQI,UE_ACK_codeblocks,UE_C,...
    UE_avg_CB_size,...
    cell_biterrors_coded,cell_biterrors_uncoded,cell_blocksize_coded,cell_blocksize_uncoded,cell_throughput_coded,cell_throughput_uncoded,...
    cell_throughput_useful,cell_FER_coded,cell_FER_uncoded,cell_used_codewords,cell_channel_error,nUE,tmp_SINR_SC_dB,tmp_PE_signal_power,tmp_PE_noise_power,tmp_Signal_plus_noise_power,tmp_Noise_power)
% Process temporary results for parallel simulations.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

for uu = 1:nUE
    results.UE_specific(uu).ACK                = UE_ACK(:,:,uu);
    results.UE_specific(uu).ACK_codeblocks     = UE_ACK_codeblocks(:,:,uu);
    results.UE_specific(uu).C                  = UE_C(:,:,uu);
    results.UE_specific(uu).avg_CB_size        = UE_avg_CB_size(:,:,uu);
    results.UE_specific(uu).rv_idx             = UE_rv_idx(:,:,uu);
    results.UE_specific(uu).RBs_assigned       = UE_RBs_assigned(:,uu);
    results.UE_specific(uu).biterrors_coded    = UE_biterrors_coded(:,:,uu);
    results.UE_specific(uu).biterrors_uncoded  = UE_biterrors_uncoded(:,:,uu);
    results.UE_specific(uu).blocksize_coded    = UE_blocksize_coded(:,:,uu);
    results.UE_specific(uu).blocksize_uncoded  = UE_blocksize_uncoded(:,:,uu);
    results.UE_specific(uu).throughput_coded   = UE_throughput_coded(:,:,uu);
    results.UE_specific(uu).throughput_uncoded = UE_throughput_uncoded(:,:,uu);
    results.UE_specific(uu).throughput_useful  = UE_throughput_useful(:,:,uu);
    results.UE_specific(uu).FER_coded          = UE_FER_coded(:,:,uu);
    results.UE_specific(uu).FER_uncoded        = UE_FER_uncoded(:,:,uu);
    results.UE_specific(uu).used_codewords     = UE_used_codewords(:,:,uu);
    results.UE_specific(uu).used_CQI           = UE_used_CQI(:,:,uu);
end

    results.cell_specific.biterrors_coded    = cell_biterrors_coded(:,:,:);
    results.cell_specific.biterrors_uncoded  = cell_biterrors_uncoded(:,:,:);
    results.cell_specific.blocksize_coded    = cell_blocksize_coded(:,:,:);
    results.cell_specific.blocksize_uncoded  = cell_blocksize_uncoded(:,:,:);
    results.cell_specific.throughput_coded   = cell_throughput_coded(:,:,:);
    results.cell_specific.throughput_uncoded = cell_throughput_uncoded(:,:,:);
    results.cell_specific.throughput_useful  = cell_throughput_useful(:,:,:);
    results.cell_specific.FER_coded          = cell_FER_coded(:,:,:);
    results.cell_specific.FER_uncoded        = cell_FER_uncoded(:,:,:);
    results.cell_specific.used_codewords     = cell_used_codewords(:,:,:);
    results.cell_specific.channel_error      = cell_channel_error(:,:,:);
    results.cell_specific.SINR_SC_dB         = tmp_SINR_SC_dB;
    results.cell_specific.PE_signal_power    = tmp_PE_signal_power;
    results.cell_specific.PE_noise_power    = tmp_PE_noise_power;
    results.cell_specific.Signal_plus_noise_power    = tmp_Signal_plus_noise_power;
    results.cell_specific.Noise_power    = tmp_Noise_power;
