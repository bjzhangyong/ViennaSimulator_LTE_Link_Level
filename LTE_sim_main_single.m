% LTE system simulator main simulator file. Check the LTE_sim_batch files
% to check how to launch the simulator.
% [] = LTE_sim_main()
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% date of creation: 2008/08/11
% last changes: 2008/09/02  Colom Ikuno HARQ is implemented
%               2008/09/04  Bosanska    changed to user_biterror and cell_biterrors... 
%                                       for multi-user implementation
%               2008/09/08  Bosanska    added structure Results.UE_specific/CELL_specific
%               2008/09/10  Bosanska    HARQ adapted for multiple users
%               2008/09/15  Bosanska    changed structure of BS_output
%               2008/09/18  Bosanska    changed the name of chan_output -> ChanMod_output
%                                       added function [ChanMod_output] = LTE_channel_matrix(ChanMod);
%               2008/10/02  Bosanska    changed structure to multiple call of function LTE_RX 
%                                       for multi-user scenario (parallel receivers)
%                                       changed function LTE_RX, new inputs -> [1x1]struct BS_UE_specific
%                                                                              [1x1]double uu
%               2008/10/21  Bosanska    changed function LTE_channel_matrix, new inputs -> [1xnUE]struct UE_output (also output) 
%                                                                                          [1x1]double uu, [1x1]double SNR
%                                       prepared loop over number of users for LTE_channel_matrix 
%                                       (multi-user scenario with different channels for users) 
%               2008/12/04  Bosanska    commented loop over number of users for LTE_channel_matrix
%               2009/02/02  Simko       scheduler plot
%               2009/03/03  Simko       multi user channel- every user has different channel
%               2009/03/09  Jcolom      Added uplink delay
%               2009/16/09  Jcolom      Splitted tracing and plotting from the main simulation loop (better code readability)
%               2009/05/15  Jcolom      Public release of the simulator (r400).
%
%
%
% By using this simulator, you agree to the license terms stated in the license agreement included with this work.
% If you are using the simulator for your scientific work, please reference:
%
% BibTeX:
% @InProceedings{EUSIPCO2009,
%   author =        {Christian Mehlf\"uhrer and Martin Wrulich and Josep Colom Ikuno and Dagmar Bosanska and Markus Rupp},
%   title =         {Simulating the Long Term Evolution Physical Layer},
%   booktitle =     {Proc. of the 17th European Signal Processing Conference (EUSIPCO 2009)},
%   month =         aug,
%   year =          2009,
%   address =       {Glasgow, Scotland},
%   note =          {accepted for publication},
% }
% 
% ASCII
% C. Mehlführer, M. Wrulich, J. C. Ikuno, D. Bosanska and M. Rupp, "Simulating the Long Term Evolution Physical Layer,"
% in Proc. of the 17th European Signal Processing Conference (EUSIPCO 2009), Aug. 2008, Glasgow, Scotland

maxStreams = 2;
simulation_results = results.simulationResults(...
    1,...
    LTE_params.nUE,...
    sum(N_subframes),...
    SNR_vec,...
    maxStreams,...
    max([UE.nRX]),...
    LTE_params.BS_config.nTx,...
    LTE_params.trace_subcarrier_SNR,...
    LTE_params.Ntot);

%% SNR and frame loops
StartTime = clock;

for SNR_i=1:size(SNR_vec,2)
    
    %% Reset of the random generators
    if LTE_params.random_channel_param_seeding
        reset(LTE_params.channel_param_RandStream,LTE_params.channel_param_seed);
    end
    if LTE_params.random_noise_seeding
        reset(LTE_params.noise_RandStream,LTE_params.noise_seed);
    end
    if LTE_params.random_data_seeding
        reset(LTE_params.data_RandStream,LTE_params.data_seed);
    end
    
    % Initialize variables that will be reused, such as the BS_output
    for bb = 1:LTE_params.nBS
        BS_output(bb) = outputs.bsOutput(LTE_params.nUE,LTE_params.Nrb,2); % Hard-coded maximum of 2 streams (obvious for LTE)
        BS_output(bb).cell_genie.SINR = zeros(LTE_params.nUE,LTE_params.Ntot);
    end
    subframe_i = 1;
    delay_counter = 0;
    
    delay = true;
    
    % Set network clock. It will tell every network element in which TTI we are right now.
    network_clock = network_elements.clock(LTE_params.FrameDur/10);
    
    % Attach the clock to each network element
    for u_=1:LTE_params.nUE*LTE_params.nBS
        UE(u_).clock = network_clock;
    end
    for b_=1:LTE_params.nBS
        BS(b_).clock = network_clock;
    end
    
    % Some further initialization
    SNR = reshape(SNR_vec(:,SNR_i,:),size(SNR_vec,1),size(SNR_vec,3));
    network_clock.reset;
    for bb = 1:LTE_params.nBS
        BS(bb).reset_HARQ_process_index;
    end
    ChanMod_output = cell(LTE_params.nBS,LTE_params.nUE*LTE_params.nBS);
    UE_input = cell(1,LTE_params.nUE*LTE_params.nBS);
    
    % Reset number of channel realizations
    for uu = 1:LTE_params.nUE*LTE_params.nBS
        UE(uu).realization_num = 0;
    end
    % Initialize uplink channel
    uplinkChannel = channels.uplinkChannel(LTE_params.uplink_delay,LTE_params.nUE);
    
    if DEBUG_LEVEL > 0
        disp('');
        disp(['*************** SNR = ' num2str(SNR(1)) 'dB, value ' num2str(SNR_i) ' of ' num2str(size(SNR_vec,2)) ' ***************']);
    end
    
    % Get channel estimation MSE
    for subs_i = 1:length(N_subframes)
        for uu=1:LTE_params.nUE*LTE_params.nBS
            UE(uu).MSE = LTE_channelestimator_MSE(10^(-SNR(uu,subs_i)/10),UE(uu).channel_autocorrelation_matrix,LTE_params,BS(1).nAtPort,UE(uu).user_speed);
        end
        while subframe_i <= sum(N_subframes(1:subs_i))
            ChanMod_output = cell(LTE_params.nBS,LTE_params.nUE*LTE_params.nBS);
            if mod(subframe_i,50)==1 && ~delay
                if DEBUG_LEVEL > 0
                    disp(['   processing subframe #' num2str(subframe_i) ' of ' num2str(sum(N_subframes))])
                    disp(['---> remaining simulation time: ' num2str(etime(clock,StartTime)/((SNR_i-1)*sum(N_subframes)+subframe_i-1)*((size(SNR_vec,2)-SNR_i)*sum(N_subframes)+sum(N_subframes)-subframe_i)/60,'%5.3f') 'min']);
                    pause(0.05);
                end
            end
            
            if ~mod(subframe_i,LTE_params.N_seed_reset) && ~delay && ~strcmp(LTE_params.ChanMod_config.time_correlation,'independent') && LTE_params.use_seed_reset
                LTE_params.channel_param_RandStream = RandStream('mt19937ar','Seed',ceil(subframe_i*13/14));
                if ~strcmp(LTE_params.ChanMod_config.type,'flat Rayleigh') && ~strcmp(LTE_params.ChanMod_config.type,'AWGN')
                    for uu=1:LTE_params.nUE*LTE_params.nBS
                        number_of_taps = ChanMod.interpolator.num_faders;
                        UE(uu).channel_coef_rosa_zheng.theta = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps)*2 -1) * pi;
                        UE(uu).channel_coef_rosa_zheng.phi = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;
                        UE(uu).channel_coef_rosa_zheng.psi = (rand(LTE_params.channel_param_RandStream,UE(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;
                    end
                end
                delay_counter = 0;
                delay = true;
                for uu = 1:LTE_params.nUE*LTE_params.nBS   % previous channels are not usable anymore for channel prediction after seed reset
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
            if LTE_params.nBS>1
                for uu = LTE_params.nUE+1:LTE_params.nUE*LTE_params.nBS
                    UE_output(uu).ACK    = true(1,2);
                    UE_output(uu).rv_idx = zeros(1,2);
                end
            end
            
            % ACK of the previous frame. If this is the first frame, set the
            % ACK to correct so that the HARQ handling generates new data
            if subframe_i==1 %|| (mod(subframe_i,LTE_params.N_seed_reset) <= LTE_params.uplink_delay && strcmp(LTE_params.ChanMod_config.time_correlation,'correlated') && LTE_params.use_seed_reset) % (number of max HARQ processes that will be used)
                for uu = 1:LTE_params.nUE
                    UE_output(uu).ACK    = true(1,2);
                    UE_output(uu).rv_idx = zeros(1,2);
                end
            end
            
            % Update current HARQ process index
            BS(1).update_current_HARQ_process(UE_output(LTE_params.connection_table(1,:)));
            
            % Generation of the channel matrix, for scheduler purposes, zero delay from RX
            switch ChanMod.type
                case 'winner_II'
                    for uu = 1:LTE_params.nUE*LTE_params.nBS % parallel channels for multi-user
                        for bb = 1:LTE_params.nBS
                            % NOTE: every user should have a different MIMO channel matrix
                            switch ChanMod.filtering
                                case 'BlockFading'
                                    [ChanMod_output{bb,uu} BS_output(bb).cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params, ChanMod, SNR(uu,subs_i), UE_output, UE(uu), uu, LTE_params.TxSymbols, subframe_i,channel{uu}(:,:,:,subframe_i),out);
                                case 'FastFading'
                                    switch ChanMod.time_correlation
                                        case 'correlated'
                                            [channel, delays, out] = LTE_winner_channel_model(LTE_params.TxSymbols,LTE_params.Arrays,out);
                                        case 'independent'
                                            [channel, delays, out] = LTE_winner_channel_model(LTE_params.TxSymbols,LTE_params.Arrays);
                                    end
                                    [ChanMod_output{bb,uu} BS_output(bb).cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params, ChanMod, SNR(uu,subs_i), UE_output, UE(uu), uu, LTE_params.TxSymbols, subframe_i,channel{uu},out);
                            end
                        end
                    end
                otherwise
                    for uu = 1:LTE_params.nUE*LTE_params.nBS % parallel channels for multi-user
                        for bb = 1:LTE_params.nBS
                            % NOTE: every user should have a different MIMO channel matrix
                            [ChanMod_output{bb,uu} BS_output(bb).cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params, ChanMod, SNR(uu,subs_i), UE_output, UE(uu), uu, LTE_params.TxSymbols, subframe_i);
                        end
                    end
            end
            
            if ~delay  % to have useful feedback values (CQI,PMI,RI) after seed reset
                %         if LTE_params.UE_config.mode == 4
                if ~LTE_params.uplink_delay
                    for uu = 1:LTE_params.nUE*LTE_params.nBS
                        [RI_tmp,PMI_tmp,CQI_tmp,CQI_bar]=LTE_feedback(BS(1).nAtPort,10^-(SNR(uu,subs_i)/10),LTE_params,ChanMod_output{uu}.genie.H_fft,UE(uu),uu,LTE_params.UE_config.mode,cqi_i);
                        UE_output(uu).RI = RI_tmp;
                        UE_output(uu).PMI = PMI_tmp-1;
                        UE_output(uu).CQI_bar = CQI_bar;
                        if ~isempty(CQI_tmp) % set the cqi value if there is no feedback
                            UE_output(uu).CQI = CQI_tmp;
                        else
                            UE_output(uu).CQI = cqi_i*ones(LTE_params.Nrb,2,LTE_params.scheduler.nLayers(uu));
                        end
                    end
                end
                if delay_counter == 0 && LTE_params.uplink_delay ~= 0
                    for uu = 1:LTE_params.nUE*LTE_params.nBS
                        UE_output(uu).CQI = ones(LTE_params.Nrb,2,LTE_params.scheduler.nLayers(uu)); % set some initial CQI values for mode 3
                    end
                end

                if LTE_params.UE_config.mode == 6 % Interference Alignment
                    P = [1 1 1];
                    sigma_n2 = 10^(-SNR_vec(SNR_i)/10);
                    [V, U, IA_channel_error(SNR_i,subframe_i)] = LTE_interference_alignment(LTE_params.IA_type, LTE_params.scheduler.nLayers, P, LTE_params.IA_thresh, LTE_params.IA_max_iterations, LTE_params, ChanMod_output, BS, UE, sigma_n2);
                end

                % Generation of the transmit signal
                % NOTE: some day there should also be a loop over basestations
                for bb = 1:LTE_params.nBS
                    if LTE_params.UE_config.mode == 6 % Interference Alignment
                        BS_output(bb).UE_signaling.MCS_and_scheduling.V = reshape(V(bb,:,:),LTE_params.Ntot,LTE_params.Nsub);
                        BS_output(bb).UE_signaling.MCS_and_scheduling.U = reshape(U(bb,:,:),LTE_params.Ntot,LTE_params.Nsub);
                    end
                    LTE_TX(LTE_params,BS(bb), UE(LTE_params.connection_table(bb,:)), BS(bb).AtPort, subframe_i,UE_output(LTE_params.connection_table(bb,:)),BS_output(bb)); % ref and sync is repeating on slot basis                    
                end
                % Convolution part
                [ChanMod_output UE_input] = LTE_channel_model(LTE_params, ChanMod, ChanMod_output, BS_output, SNR(:,subs_i), UE_input, subframe_i);
                
                % Signal receive, demod, decode...
                for uu = 1:LTE_params.nUE % parallel receivers for multi-user
                    for stream_i = 1:BS_output(1).UE_signaling(uu).MCS_and_scheduling.nCodewords % NOTE: this should be the number of streams assigned to this user by the scheduler
                        UE_output(uu).rx_data_bits{stream_i} = [];
                        UE_output(uu).rx_coded_bits{stream_i} = [];
                    end
                    UE_output(uu).PE_noise_power_subframe = nan(LTE_params.Ntot,LTE_params.Nsub);
                    UE_output(uu).PE_signal_power_subframe = nan(LTE_params.Ntot,LTE_params.Nsub);
%                     UE_output(uu).Signal_plus_noise_power = nan(UE(uu).nRX,1);
%                     UE_output(uu).Noise_power = nan(UE(uu).nRX,2);
                    % Execute receiver
                    %feedback = LTE_mini_RX(LTE_params, UE_input{uu}, SNR(uu,subs_i), subframe_i, BS_output(1), UE(1));
                    LTE_RX(LTE_params, ChanMod_output{1,uu}, UE_input{uu}, ChanMod, SNR(uu,subs_i), BS(1).AtPort, subframe_i, BS(1), UE(LTE_params.connection_table(1,:)), UE_output(LTE_params.connection_table(1,:)), BS_output(1), uu);
                end
                % SINR PostProcessing
                %PE_SINR_subframe(isnan(PE_SINR_subframe)) = 0;
                %nr_of_data_values_sub_temp = PE_SINR_subframe==0;
                %nr_of_data_values_sub = sum(nr_of_data_values_sub_temp,2);
                %PE_SINR_subcarriers = sum(PE_SINR_subframe,2)./nr_of_data_values_sub;
                
                % Put the feedback into the uplink channel (delay)
                uplinkChannel.insert_feedback(UE_output(LTE_params.connection_table(1,:)));
                
                if ~delay
                    % Process results for this TTI
                    simulation_results.process_TTI_results(BS_output(1),UE_output(LTE_params.connection_table(1,:)),subframe_i,SNR_i);
                end
                
                % Add SINR feedback (when configuration so says)
                if LTE_params.trace_subcarrier_SNR
                    % Add subcarrier SNR to the trace
                    simulation_results.cell_specific.SINR_SC_dB(:,:,subframe_i,SNR_i) = BS_output(1).cell_genie.SINR;
                end
                if delay
                    delay_counter = delay_counter+1;
                    subframe_i = subframe_i-1;
                end
                
            else
                for uu = 1:LTE_params.nUE
                    if subframe_i ~= 1
                        alphabet = BS_output.UE_signaling(uu).MCS_and_scheduling(1).CQI_params(1).modulation_order(1);
                    else
                        alphabet = 2;
                    end
                    [RI_tmp,PMI_tmp,CQI_tmp]=LTE_feedback(BS(1).nAtPort,10^-(SNR(uu,subs_i)/10),LTE_params,ChanMod_output{uu}.genie.H_fft,UE(uu),uu,LTE_params.UE_config.mode);
                    UE_output(uu).RI = RI_tmp;
                    UE_output(uu).PMI = PMI_tmp-1;
                    if ~isempty(CQI_tmp)
                        UE_output(uu).CQI = CQI_tmp;
                    else
                        UE_output(uu).CQI = cqi_i*ones(LTE_params.Nrb,2,LTE_params.scheduler.nLayers(uu));
                    end
                end
                uplinkChannel.insert_feedback(UE_output(LTE_params.connection_table(1,:)));
                delay_counter = delay_counter+1;
                subframe_i = subframe_i-1;
            end
            
            subframe_i = subframe_i+1;
        end
    end
    if DEBUG_LEVEL > 1
        UE_to_show_BLER_in_screen = 1;
        received_ACKs  = sum(simulation_results.UE_specific(UE_to_show_BLER_in_screen).ACK(:,SNR_i,1));
        scheduled_TTIs = sum(simulation_results.UE_specific(UE_to_show_BLER_in_screen).RBs_assigned(:,SNR_i)~=0);
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

%% Calculate simulation aggregates
if(strcmp(UE(1).autocorrelation_matrix_type,'estimated') && strcmp(UE(1).channel_estimation_method,'MMSE'))
    simulation_results.calculate_sim_aggregates(UE(1).realization_num_total/(UE(1).nRX * ChanMod.nTX));
else
    simulation_results.calculate_sim_aggregates(0);
end

simulation_results.SNR_vector = SNR_vec;
% simulation_results.plot_BLER_throughput;
% LTE_sim_results_scheduler_plot(BS,BS_output,subframe_i);
%% Show plots at the end of the simulation
if LTE_params.show_plots
    LTE_sim_result_plots(simulation_results);
    LTE_sim_results_scheduler_plot(BS,BS_output,subframe_i);
end

if LTE_params.UE_config.mode == 6 % Interference Alignment
    IA_channel_error = mean(IA_channel_error,2);
end