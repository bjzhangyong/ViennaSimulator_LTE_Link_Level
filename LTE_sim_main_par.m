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

clear tmp_results;
clear filename_suffix;
clear SNR_vec2;
clear StartTime;
clear maxStreams;
clear maxi;
clear num;
% clear output_filename;
clear simulation_results;


maxStreams = 2;
simulation_results = results.simulationResults(...
    LTE_params.nBS,...
    LTE_params.nUE,...
    sum(N_subframes),...
    SNR_vec,...
    maxStreams,...
    max([UE.nRX]),...
    LTE_params.BS_config.nTx,...
    LTE_params.trace_subcarrier_SNR,...
    LTE_params.Ntot);

maxi = max([UE.nRX]);

% temporary result variables for parfor
UE_res =  struct(...
    'ACK',false(sum(N_subframes),maxStreams),...
    'rv_idx',zeros(sum(N_subframes),maxStreams,'uint8'),...
    'RBs_assigned',zeros(sum(N_subframes),'uint8'),...
    'biterrors_coded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'biterrors_uncoded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'blocksize_coded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'blocksize_uncoded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'throughput_coded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'throughput_uncoded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'throughput_useful',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'FER_coded',zeros(sum(N_subframes),maxStreams,'uint16'),...
    'FER_uncoded',zeros(sum(N_subframes),maxStreams,'uint16'),...
    'used_codewords',zeros(sum(N_subframes),maxStreams,'uint16'),...
    'ACK_codeblocks',zeros(sum(N_subframes),maxStreams,'uint16'),...
    'C',zeros(sum(N_subframes),maxStreams,'uint8'),...
    'avg_CB_size',zeros(sum(N_subframes),maxStreams));

cell_res =  struct(...
    'biterrors_coded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'biterrors_uncoded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'blocksize_coded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'blocksize_uncoded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'throughput_coded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'throughput_uncoded',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'throughput_useful',zeros(sum(N_subframes),maxStreams,'uint32'),...
    'FER_coded',zeros(sum(N_subframes),maxStreams,'uint16'),...
    'FER_uncoded',zeros(sum(N_subframes),maxStreams,'uint16'),...
    'used_codewords',zeros(sum(N_subframes),maxStreams,'uint16'),...
    'channel_error',zeros(sum(N_subframes),maxi,LTE_params.BS_config.nTx),...
    'SINR_SC_dB',zeros(LTE_params.Ntot,sum(N_subframes)),...
    'PE_signal_power',nan(N_subframes,LTE_params.Ntot),...
    'PE_noise_power',nan(N_subframes,LTE_params.Ntot),...
    'Signal_plus_noise_power',nan(N_subframes,LTE_params.UE_config.nRX),...
    'Noise_power',nan(N_subframes,LTE_params.UE_config.nRX,2));

tmp_results.UE_specific = repmat(UE_res,1,LTE_params.nUE);
tmp_results.cell_specific = cell_res;
tmp_results = repmat(tmp_results,1,size(SNR_vec,2));
clear UE_res cell_res;

%% Parallel toolbox
num = matlabpool('size');
if strcmp(LTE_params.simulation_type,'parallel') && ~num
    matlabpool('open');
    %     matlabpool 1
end
par_sched = findResource();

%% SNR and frame loops
StartTime = clock;
% NOTE: changing this with a parfor does not change any functionality, but
% messes up the timer (it goes backwards), so if you really want parallel
% functionality, change the line below to a parfor
SNR_vec2 = SNR_vec;

%% File to save the progress through the SNR loop (necessary for estimation of remaining time - works only if all workers access the same harddisk)
fid = fopen('counter.m','w');
fprintf(fid,'%c','');
fclose(fid);

parfor SNR_i=1:size(SNR_vec,2)
    
    %% reset of the random generators
    if LTE_params.random_channel_param_seeding
        reset(LTE_params.channel_param_RandStream,LTE_params.channel_param_seed);
    end
    if LTE_params.random_noise_seeding
        reset(LTE_params.noise_RandStream,LTE_params.noise_seed);
    end
    if LTE_params.random_data_seeding
        reset(LTE_params.data_RandStream,LTE_params.data_seed);
    end
    
    UE_tmp = UE;
    BS_tmp = BS;
    out_tmp = out;
    LTE_params_tmp = LTE_params;
    task = getCurrentTask;
    task_ID = get(task,'ID');
    subframe_i = 1;
    delay_counter = 0;
    tmp_channel = channel;
    %     if LTE_params_tmp.UE_config.mode == 4
    delay = true;
    %     else delay = false;
    %     end
    
    % temporary result variables for parfor
    
    tmp_UE_ACK = false(sum(N_subframes),maxStreams,LTE_params_tmp.nUE);
    tmp_UE_rv_idx = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint8');
    tmp_UE_RBs_assigned = zeros(sum(N_subframes),LTE_params_tmp.nUE,'uint8');
    tmp_UE_biterrors_coded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint32');
    tmp_UE_biterrors_uncoded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint32');
    tmp_UE_blocksize_coded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint32');
    tmp_UE_blocksize_uncoded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint32');
    tmp_UE_throughput_coded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint32');
    tmp_UE_throughput_uncoded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint32');
    tmp_UE_throughput_useful = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint32');
    tmp_UE_FER_coded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint16');
    tmp_UE_FER_uncoded = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint16');
    tmp_UE_used_codewords = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint16');
    tmp_UE_used_CQI = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint8');
    tmp_UE_ACK_codeblocks = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint16');
    tmp_UE_C = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE,'uint8');
    tmp_UE_avg_CB_size = zeros(sum(N_subframes),maxStreams,LTE_params_tmp.nUE);
    
    tmp_cell_biterrors_coded = zeros(sum(N_subframes),maxStreams,'uint32');
    tmp_cell_biterrors_uncoded = zeros(sum(N_subframes),maxStreams,'uint32');
    tmp_cell_blocksize_coded = zeros(sum(N_subframes),maxStreams,'uint32');
    tmp_cell_blocksize_uncoded = zeros(sum(N_subframes),maxStreams,'uint32');
    tmp_cell_throughput_coded = zeros(sum(N_subframes),maxStreams,'uint32');
    tmp_cell_throughput_uncoded = zeros(sum(N_subframes),maxStreams,'uint32');
    tmp_cell_throughput_useful = zeros(sum(N_subframes),maxStreams,'uint32');
    tmp_cell_FER_coded = zeros(sum(N_subframes),maxStreams,'uint16');
    tmp_cell_FER_uncoded = zeros(sum(N_subframes),maxStreams,'uint16');
    tmp_cell_used_codewords = zeros(sum(N_subframes),maxStreams,'uint16');
    tmp_cell_channel_error = zeros(sum(N_subframes),maxi,LTE_params_tmp.BS_config.nTx);
    tmp_SINR_SC_dB = zeros(LTE_params_tmp.Ntot,sum(N_subframes));
    tmp_PE_signal_power = nan(N_subframes,LTE_params_tmp.Ntot);
    tmp_PE_noise_power = nan(N_subframes,LTE_params_tmp.Ntot);
    tmp_Signal_plus_noise_power = nan(N_subframes,LTE_params_tmp.UE_config.nRX);
    tmp_Noise_power = nan(N_subframes,LTE_params_tmp.UE_config.nRX,2);
    
    % Initialize variables that will be reused, such as the BS_output
    BS_output = outputs.bsOutput(LTE_params_tmp.nUE,LTE_params_tmp.Nrb,2); % Hard-coded maximum of 2 streams (obvious for LTE)
    BS_output.cell_genie.SINR = zeros(LTE_params_tmp.nUE,LTE_params_tmp.Ntot);
    
    % Set network clock. It will tell every network element in which TTI we are right now.
    network_clock = network_elements.clock(LTE_params_tmp.FrameDur/10);
    
    % Attach the clock to each network element
    
    UE_tmp = LTE_setclock(UE_tmp,network_clock);
    BS_tmp = LTE_setclock(BS_tmp,network_clock);
    
    %     for u_=1:LTE_params.nUE
    %         UE(u_,SNR_i).clock = network_clock;
    %     end
    %     for b_=1:LTE_params.nBS
    %         BS(b_).setclock(network_clock);
    %     end
    
    % Some further initialization
    SNR = reshape(SNR_vec(:,SNR_i,:),size(SNR_vec,1),size(SNR_vec,3));
    network_clock.reset;
    BS_tmp.reset_HARQ_process_index;
    ChanMod_output = cell(LTE_params.nBS,LTE_params_tmp.nUE);
    UE_input = cell(1,LTE_params.nUE*LTE_params.nBS);
    
    %     UE(:).realization_num = 0;
    %reset number of channel realizations
    for uu = 1:LTE_params_tmp.nUE
        UE_tmp(uu).realization_num = 0;
    end
    % Initialize uplink channel
    uplinkChannel = channels.uplinkChannel(LTE_params_tmp.uplink_delay,LTE_params_tmp.nUE);
    
    if DEBUG_LEVEL > 0
        disp('');
        disp(['*************** SNR = ' num2str(SNR(1)) 'dB, value ' num2str(SNR_i) ' of ' num2str(size(SNR_vec,2)) ' ***************']);
    end
    
    % get channel estimation MSE
    for subs_i = 1:length(N_subframes)
        for uu=1:LTE_params.nUE
            %         UE(uu).MSE = LTE_channelestimator_MSE(sigma_n2,UE(uu).channel_autocorrelation_matrix,LTE_params,nAtPort,UE(uu).user_speed);
            UE_tmp(uu).MSE = LTE_channelestimator_MSE(10^(-SNR(uu,subs_i)/10),UE_tmp(uu).channel_autocorrelation_matrix,LTE_params,BS(1).nAtPort,UE_tmp(uu).user_speed);
        end
        
        while subframe_i <= sum(N_subframes(1:subs_i))
            if mod(subframe_i,50)==1 && ~delay
                if DEBUG_LEVEL > 0
                    if task_ID == 1
                        job = getCurrentJob;
                        task_nr = length(get(job,'Tasks'));
                        fid = fopen('counter.m','r+');
                        A = fscanf(fid,'%d');
                        done = numel(A);
                        fclose(fid);
                        if (task_nr > size(SNR_vec2,2) - done)
                            task_nr = size(SNR_vec2,2) - done;
                        end
                        disp(['   processing subframe #' num2str(subframe_i) ' of ' num2str(sum(N_subframes))])
                        %                     disp(['---> remaining simulation time: ' num2str(etime(clock,StartTime)/((SNR_i-1)*N_subframes+subframe_i-1)*((length(SNR_vec)/task_nr-SNR_i)*N_subframes+N_subframes-subframe_i)/60,'%5.3f') 'min']);
                        disp(['---> remaining simulation time: ' num2str(etime(clock,StartTime)/60*((sum(N_subframes)*size(SNR_vec2,2))-(subframe_i-1)*task_nr-sum(N_subframes)*done)/(sum(N_subframes)*done+(subframe_i-1)*task_nr),'%5.3f') 'min']);
                        %                     pause(0.05);
                    end
                end
            end
            
            if ~mod(subframe_i,LTE_params_tmp.N_seed_reset) && ~delay && ~strcmp(LTE_params_tmp.ChanMod_config.time_correlation,'independent') && LTE_params_tmp.use_seed_reset
                LTE_params_tmp.channel_param_RandStream = RandStream('mt19937ar','Seed',ceil(subframe_i*13/14));
                if ~strcmp(LTE_params.ChanMod_config.type,'flat Rayleigh') && ~strcmp(LTE_params.ChanMod_config.type,'AWGN')
                    for uu=1:LTE_params_tmp.nUE
                        number_of_taps = ChanMod.interpolator.num_faders;
                        UE_tmp(uu).channel_coef_rosa_zheng.theta = (rand(LTE_params_tmp.channel_param_RandStream,UE_tmp(uu).nRX,BS.nTX,number_of_taps)*2 -1) * pi;
                        UE_tmp(uu).channel_coef_rosa_zheng.phi = (rand(LTE_params_tmp.channel_param_RandStream,UE_tmp(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;
                        UE_tmp(uu).channel_coef_rosa_zheng.psi = (rand(LTE_params_tmp.channel_param_RandStream,UE_tmp(uu).nRX,BS.nTX,number_of_taps,ChanMod.sin_num)*2 -1) * pi;
                    end
                end
                %             if LTE_params_tmp.UE_config.mode == 4
                delay_counter = 0;
                delay = true;
                %             end
                for uu = 1:LTE_params_tmp.nUE   % previous channels are not usable anymore for channel prediction after seed reset
                    UE_tmp(uu).previous_channels = zeros(size(UE_tmp(uu).previous_channels));
                end
            end
            
            if delay_counter == LTE_params_tmp.uplink_delay;
                delay = false;
            end
            
            % First of all, advance the network clock
            network_clock.advance_1_TTI;
            
            % Receive feedback from the previous subframe
            
            UE_output = uplinkChannel.receive_feedback;
            
            % ACK of the previous frame. If this is the first frame, set the
            % ACK to correct so that the HARQ handling generates new data
            if subframe_i==1 %|| (mod(subframe_i,LTE_params_tmp.N_seed_reset) <= LTE_params_tmp.uplink_delay && ~strcmp(LTE_params_tmp.ChanMod_config.time_correlation,'independent') && LTE_params.use_seed_reset) % (number of max HARQ processes that will be used)
                for uu = 1:LTE_params_tmp.nUE
                    UE_output(uu).ACK    = true(1,2);
                    UE_output(uu).rv_idx = zeros(1,2);
                end
            end
            
            % Update current HARQ process index
            BS_tmp.update_current_HARQ_process(UE_output);
            
            % Generation of the channel matrix, for scheduler purposes, zero delay from RX
            switch ChanMod.type
                case 'winner_II'
                    for uu = 1:LTE_params_tmp.nUE % parallel channels for multi-user
                        % NOTE: every user should have a different MIMO channel matrix
                        switch ChanMod.filtering
                            case 'BlockFading'
                                [ChanMod_output{uu} BS_output.cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params_tmp, ChanMod, SNR(uu,subs_i), UE_output, UE_tmp(uu), uu, LTE_params_tmp.TxSymbols, subframe_i,tmp_channel{uu}(:,:,:,subframe_i),out_tmp);
                            case 'FastFading'
                                switch ChanMod.time_correlation
                                    case 'correlated'
                                        [tmp_channel, delays, out_tmp] = LTE_winner_channel_model(LTE_params_tmp.TxSymbols,LTE_params_tmp.Arrays,out_tmp);
                                    case 'independent'
                                        [tmp_channel, delays, out_tmp] = LTE_winner_channel_model(LTE_params_tmp.TxSymbols,LTE_params_tmp.Arrays);
                                end
                                [ChanMod_output{uu} BS_output.cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params_tmp, ChanMod, SNR(uu,subs_i), UE_output, UE_tmp(uu), uu, LTE_params_tmp.TxSymbols, subframe_i,tmp_channel{uu},out_tmp);
                        end
                    end
                otherwise
                    for uu = 1:LTE_params_tmp.nUE % parallel channels for multi-user
                        % NOTE: every user should have a different MIMO channel matrix
                        [ChanMod_output{uu} BS_output.cell_genie.SINR(uu,:)] = LTE_channel_matrix(LTE_params_tmp, ChanMod, SNR(uu,subs_i), UE_output, UE_tmp(uu), uu, LTE_params_tmp.TxSymbols, subframe_i);
                    end
            end
            if ~delay % necessary to have useful feedback values (PMI,RI) after seed reset
                %         if LTE_params_tmp.UE_config.mode == 4
                if ~LTE_params_tmp.uplink_delay     % works just for BlockFading
                    for uu = 1:LTE_params_tmp.nUE
                        %                     if subframe_i ~= 1
                        %                         alphabet = BS_output.UE_signaling(uu).MCS_and_scheduling(1).CQI_params(1).modulation_order(1);
                        %                     else
                        %                         alphabet = 2;
                        %                     end
                        %                     [RI_tmp,PMI_tmp]=LTE_feedback_precoding(BS(1).nAtPort,10^-(SNR/10),LTE_params_tmp,alphabet,ChanMod_output{uu}.genie.H_fft,UE_tmp(uu),uu);
                        %                     [RI_tmp,PMI_tmp,CQI_tmp]=feedback_test(BS(1).nAtPort,10^-(SNR/10),LTE_params,alphabet,ChanMod_output{uu}.genie.H_fft,UE_tmp(uu),uu);
                        [RI_tmp,PMI_tmp,CQI_tmp,CQI_bar]=LTE_feedback(BS(1).nAtPort,10^-(SNR(uu,subs_i)/10),LTE_params,ChanMod_output{uu}.genie.H_fft,UE_tmp(uu),uu,LTE_params.UE_config.mode,cqi_i);
                        UE_output(uu).RI = RI_tmp;
                        UE_output(uu).PMI = PMI_tmp-1;
                        UE_output(uu).CQI_bar = CQI_bar;
                        %                     UE_output(uu).CQI = CQI_tmp;
                        if ~isempty(CQI_tmp)
                            UE_output(uu).CQI = CQI_tmp;
                        else
                            UE_output(uu).CQI = cqi_i*ones(LTE_params.Nrb,2,LTE_params.scheduler.nCodewords(uu));
                        end
                    end
                end
                if delay_counter == 0 && LTE_params_tmp.uplink_delay ~= 0
                    for uu = 1:LTE_params.nUE
                        UE_output(uu).CQI = ones(LTE_params.Nrb,2,LTE_params.scheduler.nCodewords(uu)); % set some initial CQI values for mode 3
                    end
                end
                %         end
                % Generation of the transmit signal
                % NOTE: some day there should also be a loop over basestations
                LTE_TX(LTE_params_tmp,BS_tmp, UE_tmp, BS_tmp.AtPort, subframe_i,UE_output,BS_output); % ref and sync is repeating on slot basis
                % Convolution part
                [ChanMod_output UE_input] = LTE_channel_model(LTE_params_tmp, ChanMod, ChanMod_output, BS_output, SNR(uu,subs_i),UE_input);
                % Signal receive, demod, decode...
                for uu = 1:LTE_params_tmp.nUE % parallel receivers for multi-user
                    for stream_i = 1:BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords % NOTE: this should be the number of streams assigned to this user by the scheduler
                        UE_output(uu).rx_data_bits{stream_i} = [];
                        UE_output(uu).rx_coded_bits{stream_i} = [];
                    end
                    UE_output(uu).PE_noise_power_subframe = nan(LTE_params.Ntot,LTE_params.Nsub);
                    UE_output(uu).PE_signal_power_subframe = nan(LTE_params.Ntot,LTE_params.Nsub);
                    % Execute receiver
                    LTE_RX(LTE_params_tmp, ChanMod_output{uu}, UE_input{uu}, ChanMod, SNR(uu,subs_i), BS_tmp.AtPort, subframe_i, BS_tmp, UE_tmp, UE_output, BS_output, uu);
                    
                end
                % Put the feedback into the uplink channel (delay)
                uplinkChannel.insert_feedback(UE_output);
                % Process results for this TTI
                % simulation_results.process_TTI_results(BS_output,UE_output,subframe_i,SNR_i);
                if ~delay
                    [tmp_UE_ACK(subframe_i,:,:),tmp_UE_rv_idx(subframe_i,:,:),tmp_UE_RBs_assigned(subframe_i,:),tmp_UE_biterrors_coded(subframe_i,:,:),...
                        tmp_UE_biterrors_uncoded(subframe_i,:,:),tmp_UE_blocksize_coded(subframe_i,:,:),tmp_UE_blocksize_uncoded(subframe_i,:,:),tmp_UE_throughput_coded(subframe_i,:,:),...
                        tmp_UE_FER_coded(subframe_i,:,:),tmp_UE_throughput_uncoded(subframe_i,:,:),tmp_UE_throughput_useful(subframe_i,:,:),tmp_UE_FER_uncoded(subframe_i,:,:),tmp_UE_used_codewords(subframe_i,:,:),...
                        tmp_UE_used_CQI(subframe_i,:,:),tmp_UE_ACK_codeblocks(subframe_i,:,:),tmp_UE_C(subframe_i,:,:),tmp_UE_avg_CB_size(subframe_i,:,:),...
                        tmp_cell_biterrors_coded(subframe_i,:),tmp_cell_biterrors_uncoded(subframe_i,:),tmp_cell_blocksize_coded(subframe_i,:),tmp_cell_blocksize_uncoded(subframe_i,:),...
                        tmp_cell_throughput_coded(subframe_i,:),tmp_cell_throughput_uncoded(subframe_i,:),tmp_cell_throughput_useful(subframe_i,:),tmp_cell_FER_coded(subframe_i,:),tmp_cell_FER_uncoded(subframe_i,:),...
                        tmp_cell_used_codewords(subframe_i,:),tmp_cell_channel_error(subframe_i,:,:)...
                        ] = LTE_processTTI_results(BS_output,UE_output,subframe_i,SNR_i,LTE_params_tmp.nUE,tmp_cell_biterrors_coded(subframe_i,:),...
                        tmp_cell_biterrors_uncoded(subframe_i,:),tmp_cell_blocksize_coded(subframe_i,:),tmp_cell_blocksize_uncoded(subframe_i,:),tmp_cell_throughput_coded(subframe_i,:),...
                        tmp_cell_throughput_uncoded(subframe_i,:),tmp_cell_throughput_useful(subframe_i,:),tmp_cell_FER_coded(subframe_i,:),tmp_cell_FER_uncoded(subframe_i,:),tmp_cell_used_codewords(subframe_i,:),reshape(tmp_cell_channel_error(subframe_i,:,:),maxi,LTE_params_tmp.BS_config.nTx),maxStreams);
                end
                
                % Add SINR feedback (when configuration so says)
                if LTE_params_tmp.trace_subcarrier_SNR
                    % Add subcarrier SNR to the trace
                    % simulation_results.cell_specific.SINR_SC_dB(:,subframe_i,SNR_i) = BS_output.cell_genie.SINR;
                    tmp_SINR_SC_dB(:,subframe_i) = BS_output.cell_genie.SINR;
                end
                % Postequlazation SINR post processing
                
                PE_signal_power_temp = zeros(size(UE_output(1).PE_signal_power_subframe));
                PE_noise_power_temp = zeros(size(UE_output(1).PE_signal_power_subframe));
                for uu = 1:LTE_params_tmp.nUE
                    UE_output(uu).PE_signal_power_subframe(isnan(UE_output(uu).PE_signal_power_subframe)) = 0;
                    UE_output(uu).PE_noise_power_subframe(isnan(UE_output(uu).PE_noise_power_subframe)) = 0;
                    
                    PE_signal_power_temp = PE_signal_power_temp + UE_output(uu).PE_signal_power_subframe;
                    PE_noise_power_temp = PE_noise_power_temp + UE_output(uu).PE_noise_power_subframe;
                end
                
                tmp_PE_signal_power(subframe_i,:) = sum(PE_signal_power_temp,2);
                tmp_PE_noise_power(subframe_i,:) = sum(PE_noise_power_temp,2);
                
                % estimate SNR
                tmp_Signal_plus_noise_power(subframe_i,:) = UE_output(uu).Signal_plus_noise_power;
                tmp_Noise_power(subframe_i,:,:) = UE_output(uu).Noise_power.';
                
                
                if delay
                    delay_counter = delay_counter+1;
                    subframe_i = subframe_i-1;
                end
                
            else
                for uu = 1:LTE_params_tmp.nUE
                    if subframe_i ~= 1
                        alphabet = BS_output.UE_signaling(uu).MCS_and_scheduling(1).CQI_params(1).modulation_order(1);
                    else
                        alphabet = 2;
                    end
                    [RI_tmp,PMI_tmp,CQI_tmp]=LTE_feedback(BS(1).nAtPort,10^-(SNR(uu,subs_i)/10),LTE_params,ChanMod_output{uu}.genie.H_fft,UE_tmp(uu),uu,LTE_params.UE_config.mode);
                    UE_output(uu).RI = RI_tmp;
                    UE_output(uu).PMI = PMI_tmp-1;
                    if ~isempty(CQI_tmp)
                        UE_output(uu).CQI = CQI_tmp;
                    else
                        UE_output(uu).CQI = cqi_i*ones(LTE_params.Nrb,2,LTE_params.scheduler.nCodewords(uu));
                    end
                end
                uplinkChannel.insert_feedback(UE_output);
                delay_counter = delay_counter+1;
                subframe_i = subframe_i-1;
            end
            subframe_i = subframe_i+1;
            
        end
    end
    % Set a single output variable
    tmp_results(SNR_i) = LTE_set_results(tmp_UE_ACK,tmp_UE_rv_idx,tmp_UE_RBs_assigned,tmp_UE_biterrors_coded,tmp_UE_biterrors_uncoded,tmp_UE_blocksize_coded,tmp_UE_blocksize_uncoded,...
        tmp_UE_throughput_coded,tmp_UE_throughput_uncoded,tmp_UE_throughput_useful,tmp_UE_FER_coded,tmp_UE_FER_uncoded,tmp_UE_used_codewords,tmp_UE_used_CQI,...
        tmp_UE_ACK_codeblocks,tmp_UE_C,tmp_UE_avg_CB_size,...
        tmp_cell_biterrors_coded,tmp_cell_biterrors_uncoded,...
        tmp_cell_blocksize_coded,tmp_cell_blocksize_uncoded,tmp_cell_throughput_coded,tmp_cell_throughput_uncoded,tmp_cell_throughput_useful,tmp_cell_FER_coded,tmp_cell_FER_uncoded,tmp_cell_used_codewords,tmp_cell_channel_error,LTE_params_tmp.nUE,tmp_SINR_SC_dB,tmp_PE_signal_power,tmp_PE_noise_power,tmp_Signal_plus_noise_power,tmp_Noise_power);
    
    fid = fopen('counter.m','a');
    fprintf(fid,'%d\n',SNR_i);
    fclose(fid);
end

simulation_results.SNR_vector = SNR_vec2;
simulation_results.set_TTI_results(tmp_results);
clear tmp_results;
clear U_temp;
clear BS_temp;
clear LTE_params_tmp;
clear out_tmp;

%% Calculate simulation aggregates
if(strcmp(UE(1).autocorrelation_matrix_type,'estimated') && strcmp(UE(1).channel_estimation_method,'MMSE'))
    simulation_results.calculate_sim_aggregates(UE(1).realization_num_total/(UE(1).nRX * ChanMod.nTX));
else
    simulation_results.calculate_sim_aggregates(0);
end


%% Show plots at the end of the simulation
if LTE_params.show_plots
    LTE_sim_result_plots(simulation_results);
    LTE_sim_results_scheduler_plot(BS,BS_output,sum(N_subframes));
end
if strcmp(LTE_params.simulation_type,'parallel')
    matlabpool close;
end

% To get rid of the warning when saving
clear par_sched

