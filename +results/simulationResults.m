classdef simulationResults < handle
    % Class that stores all of the simulation results. Both Cell specific and
    % UE specific.
    % Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
    % (c) 2009 by INTHFT
    % www.nt.tuwien.ac.at

    properties
        cell_specific  % Cell specific traces
        UE_specific    % UE-specific traces
        nUE
        maxStreams
        SNR_vector     % SNR vector used for this simulation
    end

    methods
        % Class constructor. Preallocate all necessary space
        function obj = simulationResults(eNodeB_count,UE_count,N_subframes,SNR_vector,maxStreams,nRx,nTx,trace_subcarrier_SNR,Ntot)

            obj.nUE = UE_count;
            obj.maxStreams = maxStreams;

            SNR_vector_length = size(SNR_vector,2);

            % Preallocate cell specific traces
            obj.cell_specific = results.cellSpecificTraces(1,1,1,1,1,Ntot,UE_count,false);
            for b_ = 1:eNodeB_count
                obj.cell_specific(b_) = results.cellSpecificTraces(N_subframes,SNR_vector_length,maxStreams,nRx,nTx,Ntot,UE_count,trace_subcarrier_SNR);
            end

            % Preallocate UE-specific traces
            obj.UE_specific = results.ueSpecificTraces(1,1,1);
            for u_ = 1:UE_count
                obj.UE_specific(u_) = results.ueSpecificTraces(N_subframes,SNR_vector_length,maxStreams);
            end
        end

        % Process results from this TTI. We will assume that only 1 eNodeB is in the simulation resutls file
        function process_TTI_results(obj,BS_output,UE_output,subframe_i,SNR_i)
            % Loop over all UEs and streams
            for uu = 1:obj.nUE
                for stream_i = 1:BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords
                    
                    if BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs
                        % Update ACK (UE traces) and biterrors (cell traces)
                        obj.UE_specific(uu).ACK(subframe_i,SNR_i,stream_i)            = UE_output(uu).ACK(stream_i);
                        obj.UE_specific(uu).ACK_codeblocks(subframe_i,SNR_i,stream_i) = UE_output(uu).ACK_codeblocks(stream_i);
                        obj.UE_specific(uu).C(subframe_i,SNR_i,stream_i)              = UE_output(uu).C(stream_i);
                        obj.UE_specific(uu).avg_CB_size(subframe_i,SNR_i,stream_i)    = mean(BS_output.UE_signaling(uu).TB_segmentation(stream_i).CB_sizes);

                        obj.UE_specific(uu).rv_idx(subframe_i,SNR_i,stream_i)         = UE_output(uu).rv_idx(stream_i);
                        obj.UE_specific(uu).RBs_assigned(subframe_i,SNR_i)            = BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs;
                        obj.UE_specific(uu).used_CQI(subframe_i,SNR_i,stream_i)       = BS_output.UE_signaling(uu).MCS_and_scheduling.cqi(stream_i);
                    else
                         % Update ACK (UE traces) and biterrors (cell traces)
                        obj.UE_specific(uu).ACK(subframe_i,SNR_i,stream_i)             = UE_output(uu).ACK(stream_i);
                        if subframe_i ~= 1
                            obj.UE_specific(uu).rv_idx(subframe_i,SNR_i,stream_i)          = obj.UE_specific(uu).rv_idx(subframe_i-1,SNR_i,stream_i);
                        else
                            obj.UE_specific(uu).rv_idx(subframe_i,SNR_i,stream_i)          = 0;
                        end
                        obj.UE_specific(uu).RBs_assigned(subframe_i,SNR_i)             = BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs;
                    end
                    
                    % Update the stats that are only meaningful when the UE has been scheduled
%                     if UE_output(uu).UE_scheduled
                    if BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs
                        % Update biterrors (UE traces)
                        if ~isempty(BS_output.genie(uu).data_bits{stream_i}) 
                            obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,stream_i)   = sum(abs(UE_output(uu).rx_data_bits{stream_i}  - BS_output.genie(uu).data_bits{stream_i}));
                        else
                            obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,stream_i)   = 0;
                        end
                        obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,stream_i) = sum(abs(UE_output(uu).rx_coded_bits{stream_i} - BS_output.genie(uu).sent_bits{stream_i}));
                        
                        % Update blocksize (UE traces)
                        obj.UE_specific(uu).blocksize_coded(subframe_i,SNR_i,stream_i)   = length(UE_output(uu).rx_data_bits{stream_i});
                        obj.UE_specific(uu).blocksize_uncoded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_coded_bits{stream_i});
                        
                        % Update FER and throughput (UE traces)
                        
                        % Coded
                        if UE_output(uu).ACK(stream_i)
                            obj.UE_specific(uu).throughput_coded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_data_bits{stream_i});
                            obj.UE_specific(uu).throughput_useful(subframe_i,SNR_i,stream_i) = BS_output.UE_signaling(uu).MCS_and_scheduling.N_used_bits(stream_i);
                            obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,stream_i) = 0;
                        else
                            obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,stream_i) = 1;
                        end
                        
                        % Uncoded
                        if (obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,stream_i)==0)
                            obj.UE_specific(uu).throughput_uncoded(subframe_i,SNR_i,stream_i) = length(UE_output(uu).rx_coded_bits{stream_i});
                            obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,stream_i)        = 0;
                        else
                            obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,stream_i) = 1;
                        end
                        
                        % Update what codewords were used (valid positions in the traces, 0 if no RBs were allocated)
                        obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,stream_i) = 1;
                    else
                        % Update what codewords were used (valid positions in the traces, 0 if no RBs were allocated)
                        obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,stream_i) = 0;
                    end
                end

                % Update cell coded and uncoded bit errors
                obj.cell_specific.biterrors_coded(subframe_i,SNR_i,:)   = obj.cell_specific.biterrors_coded(subframe_i,SNR_i,:)   + obj.UE_specific(uu).biterrors_coded(subframe_i,SNR_i,:);
                obj.cell_specific.biterrors_uncoded(subframe_i,SNR_i,:) = obj.cell_specific.biterrors_uncoded(subframe_i,SNR_i,:) + obj.UE_specific(uu).biterrors_uncoded(subframe_i,SNR_i,:);

                % Update blocksize (cell traces)
                obj.cell_specific.blocksize_coded(subframe_i,SNR_i,:)   = obj.cell_specific.blocksize_coded(subframe_i,SNR_i,:)   + obj.UE_specific(uu).blocksize_coded(subframe_i,SNR_i,:);
                obj.cell_specific.blocksize_uncoded(subframe_i,SNR_i,:) = obj.cell_specific.blocksize_uncoded(subframe_i,SNR_i,:) + obj.UE_specific(uu).blocksize_uncoded(subframe_i,SNR_i,:);

                % Update FER and throughput (cell traces)
                obj.cell_specific.throughput_coded(subframe_i,SNR_i,:)  = obj.cell_specific.throughput_coded(subframe_i,SNR_i,:) + obj.UE_specific(uu).throughput_coded(subframe_i,SNR_i,:);
                obj.cell_specific.throughput_uncoded(subframe_i,SNR_i,:)= obj.cell_specific.throughput_uncoded(subframe_i,SNR_i,:) + obj.UE_specific(uu).throughput_uncoded(subframe_i,SNR_i,:);
                obj.cell_specific.throughput_useful(subframe_i,SNR_i,:) = obj.cell_specific.throughput_useful(subframe_i,SNR_i,:) + obj.UE_specific(uu).throughput_useful(subframe_i,SNR_i,:);
                obj.cell_specific.FER_coded(subframe_i,SNR_i,:)         = obj.cell_specific.FER_coded(subframe_i,SNR_i,:)   + uint16(obj.UE_specific(uu).FER_coded(subframe_i,SNR_i,:));
                obj.cell_specific.FER_uncoded(subframe_i,SNR_i,:)       = obj.cell_specific.FER_uncoded(subframe_i,SNR_i,:) + uint16(obj.UE_specific(uu).FER_uncoded(subframe_i,SNR_i,:));
                
                % The number of codewords is the maximum (bitwise OR) of the used codewords matrix
                obj.cell_specific.used_codewords(subframe_i,SNR_i,:) = obj.cell_specific.used_codewords(subframe_i,SNR_i,:) + uint16(obj.UE_specific(uu).used_codewords(subframe_i,SNR_i,:));
                
                % If 0 RBs were allocated => no transmission  -> no channel error calculation for this user
                if(sum(BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs))
                    channel_error_help = zeros(size(obj.cell_specific.channel_error(subframe_i,SNR_i,:,:)));
                    channel_error_help(1,1,:,:) = UE_output(uu).channel_estimation_error;
                    obj.cell_specific.channel_error(subframe_i,SNR_i,:,:) = obj.cell_specific.channel_error(subframe_i,SNR_i,:,:) + channel_error_help;
                end
                
                % The estimation error of the carrier frequency offset
                if ~isempty(UE_output(uu).freq_offset_est)
                    obj.UE_specific(uu).freq_offset_est_error(subframe_i,SNR_i) = UE_output(uu).freq_offset_est.error;
                    obj.UE_specific(uu).freq_offset_est_frac(subframe_i,SNR_i) = UE_output(uu).freq_offset_est.frac;
                    obj.UE_specific(uu).freq_offset_est_int(subframe_i,SNR_i) = UE_output(uu).freq_offset_est.int;
                    obj.UE_specific(uu).freq_offset_est_res(subframe_i,SNR_i) = UE_output(uu).freq_offset_est.res;
                end
            end
            obj.cell_specific.channel_error(subframe_i,SNR_i,:,:) = obj.cell_specific.channel_error(subframe_i,SNR_i,:,:)/obj.nUE;
            
            PE_signal_power_temp = zeros(size(UE_output(uu).PE_signal_power_subframe));
            PE_noise_power_temp = zeros(size(UE_output(uu).PE_signal_power_subframe));
            for uu = 1:obj.nUE
                UE_output(uu).PE_signal_power_subframe(isnan(UE_output(uu).PE_signal_power_subframe)) = 0;
                UE_output(uu).PE_noise_power_subframe(isnan(UE_output(uu).PE_noise_power_subframe)) = 0;
                
                PE_signal_power_temp = PE_signal_power_temp + UE_output(uu).PE_signal_power_subframe;
                PE_noise_power_temp = PE_noise_power_temp + UE_output(uu).PE_noise_power_subframe;
            end
          
            obj.cell_specific.PE_signal_power(subframe_i,SNR_i,:)           = sum(PE_signal_power_temp,2);
            obj.cell_specific.PE_noise_power(subframe_i,SNR_i,:)            = sum(PE_noise_power_temp,2);
            if ~isempty(UE_output(uu).Signal_plus_noise_power)
                obj.cell_specific.Signal_plus_noise_power(subframe_i,SNR_i,:)   = UE_output(uu).Signal_plus_noise_power;
            end
            if ~isempty(UE_output(uu).Noise_power)
                obj.cell_specific.Noise_power(subframe_i,SNR_i,:,:)             = UE_output(uu).Noise_power.';
            end
            
        end

        % Calculate simulations aggregates.
        function calculate_sim_aggregates(obj,elements_to_remove)
            %remove results, which havent used estimated autocorrelations matrix
            obj.cell_specific.biterrors_coded(1:elements_to_remove,:,:)             = [];
            obj.cell_specific.biterrors_uncoded(1:elements_to_remove,:,:)           = [];
            obj.cell_specific.FER_coded(1:elements_to_remove,:,:)                   = [];
            obj.cell_specific.channel_error(1:elements_to_remove,:,:,:)             = [];
            obj.cell_specific.blocksize_coded(1:elements_to_remove,:,:)             = [];
            obj.cell_specific.blocksize_uncoded(1:elements_to_remove,:,:)           = [];
            obj.cell_specific.used_codewords(1:elements_to_remove,:,:)              = [];
            obj.cell_specific.throughput_coded(1:elements_to_remove,:,:)            = [];
            obj.cell_specific.throughput_uncoded(1:elements_to_remove,:,:)          = [];
            obj.cell_specific.FER_uncoded(1:elements_to_remove,:,:)                 = [];
            obj.cell_specific.PE_signal_power(1:elements_to_remove,:,:)             = []; 
            obj.cell_specific.PE_noise_power(1:elements_to_remove,:,:)              = [];
            obj.cell_specific.Signal_plus_noise_power(1:elements_to_remove,:,:)     = []; 
            obj.cell_specific.Noise_power(1:elements_to_remove,:,:,:)               = []; 
            
            % Cell specific results
            obj.cell_specific.BER_coded   = squeeze(sum(obj.cell_specific.biterrors_coded,1)   ./ sum(obj.cell_specific.blocksize_coded,1));
            obj.cell_specific.BER_uncoded = squeeze(sum(obj.cell_specific.biterrors_uncoded,1) ./ sum(obj.cell_specific.blocksize_uncoded,1));
            obj.cell_specific.BLER        = squeeze(sum(obj.cell_specific.FER_coded,1) ./ sum(obj.cell_specific.used_codewords,1));
            
            obj.cell_specific.BER_coded_overall   = squeeze(sum(sum(obj.cell_specific.biterrors_coded,1),3)   ./ sum(sum(obj.cell_specific.blocksize_coded,1),3));
            obj.cell_specific.BER_uncoded_overall = squeeze(sum(sum(obj.cell_specific.biterrors_uncoded,1),3) ./ sum(sum(obj.cell_specific.blocksize_uncoded,1),3));
            obj.cell_specific.BLER_overall        = squeeze(sum(sum(obj.cell_specific.FER_coded,1),3) ./ sum(sum(obj.cell_specific.used_codewords,1),3));
            obj.cell_specific.MSE_overall         = mean(mean(mean(obj.cell_specific.channel_error,4),3),1);

            obj.cell_specific.PE_SINR_overall     = 10*log10(squeeze(sum(obj.cell_specific.PE_signal_power,1)./sum(obj.cell_specific.PE_noise_power,1)));
            
            % SNR estimation
            noise_power = sum(sum(mean(obj.cell_specific.Noise_power,4),3),1);
            signal_plus_noise_power = sum(sum(obj.cell_specific.Signal_plus_noise_power,3),1);
            signal_power = signal_plus_noise_power-noise_power;
            obj.cell_specific.SNR_estimated       = 10*log10(signal_power./noise_power);
            
            % UE-specific results
            for uu = 1:obj.nUE
                obj.UE_specific(uu).biterrors_coded(1:elements_to_remove,:,:)   = [];
                obj.UE_specific(uu).biterrors_uncoded(1:elements_to_remove,:,:) = [];
                obj.UE_specific(uu).FER_coded(1:elements_to_remove,:,:)         = [];
                obj.UE_specific(uu).blocksize_coded(1:elements_to_remove,:,:)   = [];
                obj.UE_specific(uu).blocksize_uncoded(1:elements_to_remove,:,:) = [];
                obj.UE_specific(uu).used_codewords(1:elements_to_remove,:,:)    = [];
                
                obj.UE_specific(uu).BER_coded    = sum(obj.UE_specific(uu).biterrors_coded,1)   ./ sum(obj.UE_specific(uu).blocksize_coded,1);
                obj.UE_specific(uu).BER_uncoded  = sum(obj.UE_specific(uu).biterrors_uncoded,1) ./ sum(obj.UE_specific(uu).blocksize_uncoded,1);
                obj.UE_specific(uu).BLER         = squeeze(sum(obj.UE_specific(uu).FER_coded,1) ./ sum(obj.UE_specific(uu).used_codewords,1));
                
                obj.UE_specific(uu).BER_coded_overall   = squeeze(sum(sum(obj.UE_specific(uu).biterrors_coded,1),3)   ./ sum(sum(obj.UE_specific(uu).blocksize_coded,1),3));
                obj.UE_specific(uu).BER_uncoded_overall = squeeze(sum(sum(obj.UE_specific(uu).biterrors_uncoded,1),3) ./ sum(sum(obj.UE_specific(uu).blocksize_uncoded,1),3));
                obj.UE_specific(uu).BLER_overall        = squeeze(sum(sum(obj.UE_specific(uu).FER_coded,1),3) ./ sum(sum(obj.UE_specific(uu).used_codewords,1),3));
                obj.UE_specific(uu).MSE_freq_offset     = squeeze(mean(obj.UE_specific(uu).freq_offset_est_error.^2));%obj.UE_specific(uu).freq_offset_est_error;%
            end
        end
        
		% Dumps the output struct from the parallel simulation mode into the results object
        function set_TTI_results(obj,tmp_results)
            for ii = 1:size(obj.SNR_vector,2)
                for uu = 1:obj.nUE
                    obj.UE_specific(uu).ACK(:,ii,:)                = tmp_results(ii).UE_specific(uu).ACK(:,:);
                    obj.UE_specific(uu).ACK_codeblocks(:,ii,:)     = tmp_results(ii).UE_specific(uu).ACK_codeblocks(:,:);
                    obj.UE_specific(uu).C(:,ii,:)                  = tmp_results(ii).UE_specific(uu).C(:,:);
                    obj.UE_specific(uu).avg_CB_size(:,ii,:)        = tmp_results(ii).UE_specific(uu).avg_CB_size;
                    obj.UE_specific(uu).rv_idx(:,ii,:)             = tmp_results(ii).UE_specific(uu).rv_idx(:,:);
                    obj.UE_specific(uu).RBs_assigned(:,ii)         = tmp_results(ii).UE_specific(uu).RBs_assigned(:);
                    obj.UE_specific(uu).biterrors_coded(:,ii,:)    = tmp_results(ii).UE_specific(uu).biterrors_coded(:,:);
                    obj.UE_specific(uu).biterrors_uncoded(:,ii,:)  = tmp_results(ii).UE_specific(uu).biterrors_uncoded(:,:);
                    obj.UE_specific(uu).blocksize_coded(:,ii,:)    = tmp_results(ii).UE_specific(uu).blocksize_coded(:,:);
                    obj.UE_specific(uu).blocksize_uncoded(:,ii,:)  = tmp_results(ii).UE_specific(uu).blocksize_uncoded(:,:);
                    obj.UE_specific(uu).throughput_coded(:,ii,:)   = tmp_results(ii).UE_specific(uu).throughput_coded(:,:);
                    obj.UE_specific(uu).throughput_uncoded(:,ii,:) = tmp_results(ii).UE_specific(uu).throughput_uncoded(:,:);
                    obj.UE_specific(uu).throughput_useful(:,ii,:)  = tmp_results(ii).UE_specific(uu).throughput_useful(:,:);
                    obj.UE_specific(uu).FER_coded(:,ii,:)          = tmp_results(ii).UE_specific(uu).FER_coded(:,:);
                    obj.UE_specific(uu).FER_uncoded(:,ii,:)        = tmp_results(ii).UE_specific(uu).FER_uncoded(:,:);
                    obj.UE_specific(uu).used_codewords(:,ii,:)     = tmp_results(ii).UE_specific(uu).used_codewords(:,:);
                    obj.UE_specific(uu).used_CQI(:,ii,:)           = tmp_results(ii).UE_specific(uu).used_CQI(:,:);
                end
                obj.cell_specific.biterrors_coded(:,ii,:)           = tmp_results(ii).cell_specific.biterrors_coded(:,:);
                obj.cell_specific.biterrors_uncoded(:,ii,:)         = tmp_results(ii).cell_specific.biterrors_uncoded(:,:);
                obj.cell_specific.blocksize_coded(:,ii,:)           = tmp_results(ii).cell_specific.blocksize_coded(:,:);
                obj.cell_specific.blocksize_uncoded(:,ii,:)         = tmp_results(ii).cell_specific.blocksize_uncoded(:,:);
                obj.cell_specific.throughput_coded(:,ii,:)          = tmp_results(ii).cell_specific.throughput_coded(:,:);
                obj.cell_specific.throughput_uncoded(:,ii,:)        = tmp_results(ii).cell_specific.throughput_uncoded(:,:);
                obj.cell_specific.throughput_useful(:,ii,:)         = tmp_results(ii).cell_specific.throughput_useful(:,:);
                obj.cell_specific.FER_coded(:,ii,:)                 = tmp_results(ii).cell_specific.FER_coded(:,:);
                obj.cell_specific.FER_uncoded(:,ii,:)               = tmp_results(ii).cell_specific.FER_uncoded(:,:);
                obj.cell_specific.used_codewords(:,ii,:)            = tmp_results(ii).cell_specific.used_codewords(:,:);
                obj.cell_specific.channel_error(:,ii,:,:)           = tmp_results(ii).cell_specific.channel_error(:,:,:);
                obj.cell_specific.PE_signal_power(:,ii,:)           = tmp_results(ii).cell_specific.PE_signal_power(:,:);
                obj.cell_specific.PE_noise_power(:,ii,:)            = tmp_results(ii).cell_specific.PE_noise_power(:,:);
                obj.cell_specific.Signal_plus_noise_power(:,ii,:)   = tmp_results(ii).cell_specific.Signal_plus_noise_power(:,:);
                obj.cell_specific.Noise_power(:,ii,:,:)             = tmp_results(ii).cell_specific.Noise_power(:,:,:);
                
            end
        end
        
        function plot_BLER_throughput(obj,varargin)
            
            if isempty(varargin)
                first_figure = figure;
                second_figure = figure;
            else
                first_figure  = varargin{1}(1);
                second_figure = varargin{1}(2);
            end
            
            first_axes = axes('Parent',first_figure);
            second_axes = axes('Parent',second_figure);
            
            N_subframes = size(obj.cell_specific.used_codewords,1);
            plot(first_axes,obj.SNR_vector,obj.cell_specific.BLER_overall);
            plot(second_axes,obj.SNR_vector,squeeze(sum(sum(obj.cell_specific.throughput_coded,1),3))/(N_subframes*1e-3)/1e6);
            
            set(first_axes,'YScale','log');
            ylabel(first_axes,'BLER');
            xlabel(first_axes,'SNR [dB]');
            ylim(first_axes,[1e-3 1]);
            grid(first_axes,'on');
            
            
            ylabel(second_axes,'throughput [Mbps]');
            xlabel(second_axes,'SNR [dB]');
            grid(second_axes,'on');
            
        end
    end
end