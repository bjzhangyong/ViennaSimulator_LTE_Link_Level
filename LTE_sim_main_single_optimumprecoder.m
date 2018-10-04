% This file is necessary for WSA_schwarz_batch
% It chooses for every channel realization the optimal precoding matrix
% (from the LTE codebook) and uses this one for transmission

maxStreams = 2;
simulation_results = results.simulationResults(...
    LTE_params.nBS,...
    LTE_params.nUE,...
    N_subframes,...
    SNR_vec,...
    maxStreams,...
    max([UE.nRX]),...
    LTE_params.BS_config.nTx,...
    LTE_params.trace_subcarrier_SNR,...
    LTE_params.Ntot);


%% SNR and frame loops
StartTime = clock;
% NOTE: changing this with a parfor does not change any functionality, but
% messes up the timer (it goes backwards), so if you really want parallel
% functionality, change the line below to a parfor

for SNR_i=1:length(SNR_vec)
    % Initialize variables that will be reused, such as the BS_output
    BS_output = outputs.bsOutput(LTE_params.nUE,LTE_params.Nrb,2); % Hard-coded maximum of 2 streams (obvious for LTE)
    BS_output.cell_genie.SINR = zeros(LTE_params.nUE,LTE_params.Ntot);
    subframe_i = 1;
    delay_counter = 0;
      
    if LTE_params.UE_config.mode == 4
        delay = true;
    else delay = false;
    end
    % Set network clock. It will tell every network element in which TTI we are right now.
    network_clock = network_elements.clock(LTE_params.FrameDur/10);

    % Attach the clock to each network element
    for u_=1:LTE_params.nUE
        UE(u_).clock = network_clock;
    end
    for b_=1:LTE_params.nBS
        BS(b_).clock = network_clock;
    end

    % Some further initialization
    SNR = SNR_vec(SNR_i);
    network_clock.reset;
    BS.reset_HARQ_process_index;
    ChanMod_output = cell(1,LTE_params.nUE);
    %reset number of channel realizations
    for uu = 1:LTE_params.nUE
        UE(uu).realization_num = 0;
    end
    % Initialize uplink channel
    uplinkChannel = channels.uplinkChannel(LTE_params.uplink_delay,LTE_params.nUE);

    if DEBUG_LEVEL > 0
        disp('');
        disp(['*************** SNR = ' num2str(SNR) 'dB, value ' num2str(SNR_i) ' of ' num2str(length(SNR_vec)) ' ***************']);
    end

    while subframe_i <= N_subframes 
        errors = 10^10;
        if mod(subframe_i,50)==1 && ~delay
            if DEBUG_LEVEL > 0
                disp(['   processing subframe #' num2str(subframe_i) ' of ' num2str(N_subframes)])
                disp(['---> remaining simulation time: ' num2str(etime(clock,StartTime)/((SNR_i-1)*N_subframes+subframe_i-1)*((length(SNR_vec)-SNR_i)*N_subframes+N_subframes-subframe_i)/60,'%5.3f') 'min']);
                pause(0.05);
            end
        end

        if ~mod(subframe_i,LTE_params.N_seed_reset) && ~delay && ~strcmp(LTE_params.ChanMod_config.time_correlation,'independent')
            LTE_params.channel_param_RandStream = RandStream('mt19937ar','Seed',ceil(subframe_i*13/14));
            for uu=1:LTE_params.nUE
                number_of_taps = ChanMod.interpolator.num_faders;
                UE(uu).channel_coef_rosa_zhneg.theta = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps)*2 -1) * pi;
                UE(uu).channel_coef_rosa_zhneg.phi = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;
                UE(uu).channel_coef_rosa_zhneg.psi = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;
            end
            if LTE_params.UE_config.mode == 4
                delay_counter = 0;
                delay = true;
            end
            for uu = 1:LTE_params.nUE   % previous channels are not usable anymore for channel prediction after seed reset
                UE(uu).previous_channels = zeros(size(UE(uu).previous_channels));            
            end
        end

        if delay_counter == LTE_params.uplink_delay;
            delay = false;
        end

        % First of all, advance the network clock
        network_clock.advance_1_TTI;

        % Receive feedback from the previous subframe
        UE_output = uplinkChannel.receive_feedback;

        % ACK of the previous frame. If this is the first frame, set the
        % ACK to correct so that the HARQ handling generates new data
        if subframe_i==1 || (mod(subframe_i,LTE_params.N_seed_reset) <= LTE_params.uplink_delay && ~strcmp(LTE_params.ChanMod_config.time_correlation,'independent'))% (number of max HARQ processes that will be used)
            for uu = 1:LTE_params.nUE
                UE_output(uu).ACK    = true(1,2);
                UE_output(uu).rv_idx = zeros(1,2);
            end
        end

        % Update current HARQ process index
        BS.update_current_HARQ_process(UE_output);
        
        % Generation of the channel matrix, for scheduler purposes, zero delay from RX
        switch ChanMod.type
            case 'winner_II'
                for uu = 1:LTE_params.nUE % parallel channels for multi-user
                    % NOTE: every user should have a different MIMO channel matrix
                    switch ChanMod.filtering
                        case 'BlockFading'
                            [ChanMod_output{uu} BS_output.cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params, ChanMod, SNR, UE_output, UE(uu), uu, LTE_params.TxSymbols, subframe_i,channel{uu}(:,:,:,subframe_i),out);
                        case 'FastFading'
                            switch ChanMod.time_correlation
                                case 'correlated'
                                    [channel, delays, out] = LTE_winner_channel_model(LTE_params.TxSymbols,LTE_params.Arrays,out);
                                case 'independent'
                                    [channel, delays, out] = LTE_winner_channel_model(LTE_params.TxSymbols,LTE_params.Arrays);
                            end
                            [ChanMod_output{uu} BS_output.cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params, ChanMod, SNR, UE_output, UE(uu), uu, LTE_params.TxSymbols, subframe_i,channel{uu},out);
                    end
                end
            otherwise
                for uu = 1:LTE_params.nUE % parallel channels for multi-user
                    % NOTE: every user should have a different MIMO channel matrix
                    [ChanMod_output{uu} BS_output.cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params, ChanMod, SNR, UE_output, UE(uu), uu, LTE_params.TxSymbols, subframe_i);
                end
        end
        
        if ~delay % necessary to have useful feedback values (PMI,RI) after seed reset
        for ttt = 2:i_max   
        if LTE_params.UE_config.mode == 4
            if ~LTE_params.uplink_delay && UE(uu).PMI_fb
                for uu = 1:LTE_params.nUE 
                    if subframe_i ~= 1
                        alphabet = BS_output.UE_signaling(uu).MCS_and_scheduling(1).CQI_params(1).modulation_order(1);
                    else
                        alphabet = 2;
                    end    
%                     if ~isempty(BS_output.UE_signaling(uu).MCS_and_scheduling)
%                         mapping = BS_output.UE_signaling.MCS_and_scheduling.UE_mapping;
%                     else
%                         mapping = ones(LTE_params.Nrb,2);
%                     end
%                     [RI_tmp,PMI_tmp]=LTE_feedback_precoding(BS(1).nAtPort,10^-(SNR/10),LTE_params,alphabet,ChanMod_output{uu}.genie.H_fft,UE(uu));
%                     UE_output(uu).RI = RI_tmp;
%                     UE_output(uu).PMI = PMI_tmp-1;
%                         UE_output(uu).CQI_feedback = CQI;
                    if i_max == 4
                      UE_output(uu).RI = 1;
                    else
                      UE_output(uu).RI = 2;
                    end
                      UE_output(uu).PMI = (ttt-1)*ones(LTE_params.Nrb,2);
                end
            end
        end

        % Generation of the transmit signal
        % NOTE: some day there should also be a loop over basestations
        LTE_TX(LTE_params,BS, UE, BS.AtPort, subframe_i,UE_output,BS_output); % ref and sync is repeating on slot basis

        % Convolution part
        for uu = 1:LTE_params.nUE
            [ChanMod_output{uu}] = LTE_channel_model(LTE_params, ChanMod, ChanMod_output{uu}, BS_output, SNR);
        end
        
        % Signal receive, demod, decode...
        for uu = 1:LTE_params.nUE % parallel receivers for multi-user
            for stream_i = 1:BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords % NOTE: this should be the number of streams assigned to this user by the scheduler
                UE_output(uu).rx_data_bits{stream_i} = [];
                UE_output(uu).rx_coded_bits{stream_i} = [];
            end
            
            % Execute receiver
            LTE_RX(LTE_params, ChanMod_output{uu}, ChanMod, SNR, BS.AtPort, subframe_i, BS, UE, UE_output, BS_output, uu);
        end
        if BS_output.UE_signaling(1).MCS_and_scheduling.nCodewords == 2
            errors_act = mean([sum(abs(UE_output(1).rx_coded_bits{2} - BS_output.genie(1).sent_bits{2})),sum(abs(UE_output(1).rx_coded_bits{1} - BS_output.genie(1).sent_bits{1}))]);
        else
            errors_act = sum(abs(UE_output(1).rx_coded_bits{1} - BS_output.genie(1).sent_bits{1}));
        end
        if (errors > errors_act)
            opt = ttt;
            errors = errors_act;
            simulation_results.cell_specific.FER_coded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.FER_coded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.FER_uncoded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.FER_uncoded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.throughput_coded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.throughput_coded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.throughput_uncoded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.throughput_uncoded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.biterrors_coded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.biterrors_coded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.biterrors_uncoded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.biterrors_uncoded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.blocksize_coded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.blocksize_coded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.blocksize_uncoded(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.blocksize_uncoded(subframe_i,SNR_i,:);
            simulation_results.cell_specific.used_codewords(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.used_codewords(subframe_i,SNR_i,:);
            simulation_results.cell_specific.channel_error(subframe_i,SNR_i,:) = 0*simulation_results.cell_specific.channel_error(subframe_i,SNR_i,:);
            % Put the feedback into the uplink channel (delay)
            uplinkChannel.insert_feedback(UE_output);

            % Process results for this TTI
            simulation_results.process_TTI_results(BS_output,UE_output,subframe_i,SNR_i);

            % Add SINR feedback (when configuration so says)
            if LTE_params.trace_subcarrier_SNR
                % Add subcarrier SNR to the trace
                simulation_results.cell_specific.SINR_SC_dB(:,:,subframe_i,SNR_i) = BS_output.cell_genie.SINR;
            end
        end
        end
%         opt-1
        else   
            for uu = 1:LTE_params.nUE 
                    if subframe_i ~= 1
                        alphabet = BS_output.UE_signaling(uu).MCS_and_scheduling(1).CQI_params(1).modulation_order(1);
                    else
                        alphabet = 2;
                    end  
                    [RI_tmp,PMI_tmp]=LTE_feedback_precoding(BS(1).nAtPort,10^-(SNR/10),LTE_params,alphabet,ChanMod_output{uu}.genie.H_fft,UE(uu),uu);
                    UE_output(uu).RI = RI_tmp;
                    UE_output(uu).PMI = PMI_tmp-1;
            end
            uplinkChannel.insert_feedback(UE_output);
            delay_counter = delay_counter+1;
            subframe_i = subframe_i-1;
        end
        
        subframe_i = subframe_i+1;
    end
    
    if DEBUG_LEVEL > 1
        UE_to_show_BLER_in_screen = 1;
        received_ACKs  = sum(simulation_results.UE_specific(UE_to_show_BLER_in_screen).ACK(:,SNR_i,1));
        scheduled_TTIs = sum(simulation_results.UE_specific(UE_to_show_BLER_in_screen).RBs_assigned(:,SNR_i)~=0);
        disp([sprintf('   BLER UE%d, stream 1: ',UE_to_show_BLER_in_screen) num2str(1 - received_ACKs/scheduled_TTIs,'%3.2f')]);
    end
    
    % Saving and reusing the channel model trace is supported only for
    % non-parallel simulations, as using parallel Matlab instances will
    % break the method used to repeat the noise realizations
    if LTE_params.store_channel_trace
        channel_matrix_tracefile = LTE_params.channel_matrix_trace;
        if DEBUG_LEVEL > 4
            channel_matrix_tracefile.plot_trace;
        end
        save(LTE_params.channel_matrix_tracefile,'channel_matrix_tracefile');
        clear channel_matrix_tracefile;
    end
end

% if strcmp(LTE_params.simulation_type,'parallel')
%     matlabpool close;
% end

%% Calculate simulation aggregates
if(strcmp(UE(1).autocorrelation_matrix_type,'estimated') && strcmp(UE(1).channel_estimation_method,'MMSE'))
    simulation_results.calculate_sim_aggregates(UE(1).realization_num_total/(UE(1).nRX * ChanMod.nTX));
else
    simulation_results.calculate_sim_aggregates(0);
end

simulation_results.SNR_vector = SNR_vec;
% simulation_results.plot_BLER_throughput;

%% Show plots at the end of the simulation
if LTE_params.show_plots
    LTE_sim_result_plots(simulation_results);
    LTE_sim_results_scheduler_plot(BS,BS_output,subframe_i);
end
