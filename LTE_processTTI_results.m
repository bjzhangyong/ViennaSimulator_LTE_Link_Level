function [UE_ACK,UE_rv_idx,UE_RBs_assigned,UE_biterrors_coded,UE_biterrors_uncoded,UE_blocksize_coded,UE_blocksize_uncoded,UE_throughput_coded,...
    UE_FER_coded,UE_throughput_uncoded,UE_throughput_useful,UE_FER_uncoded,UE_used_codewords,UE_used_CQI,UE_ACK_codeblocks,UE_C,UE_avg_CB_size,...
    cell_biterrors_coded,cell_biterrors_uncoded,cell_blocksize_coded,cell_blocksize_uncoded,cell_throughput_coded,cell_throughput_uncoded,cell_throughput_useful,...
    cell_FER_coded,cell_FER_uncoded,cell_used_codewords,cell_channel_error...
    ] = LTE_processTTI_results(BS_output,UE_output,subframe_i,SNR_i,nUE,cell_biterrors_coded,cell_biterrors_uncoded,cell_blocksize_coded,...
    cell_blocksize_uncoded,cell_throughput_coded,cell_throughput_uncoded,cell_throughput_useful,cell_FER_coded,cell_FER_uncoded,cell_used_codewords,cell_channel_error,maxstreams)
% Process temporary results for parallel simulations.
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

% Initialization
UE_ACK                 = false(maxstreams,nUE);
UE_rv_idx              = zeros(maxstreams,nUE,'uint8');
UE_RBs_assigned        = zeros(nUE,1,'uint8');
UE_biterrors_coded     = zeros(maxstreams,nUE,'uint32');
UE_biterrors_uncoded   = zeros(maxstreams,nUE,'uint32');
UE_blocksize_coded     = zeros(maxstreams,nUE,'uint32');
UE_blocksize_uncoded   = zeros(maxstreams,nUE,'uint32');
UE_throughput_coded    = zeros(maxstreams,nUE,'uint32');
UE_throughput_uncoded  = zeros(maxstreams,nUE,'uint32');
UE_throughput_useful   = zeros(maxstreams,nUE,'uint32');
UE_FER_coded           = zeros(maxstreams,nUE,'uint16');
UE_FER_uncoded         = zeros(maxstreams,nUE,'uint16');
UE_used_codewords      = zeros(maxstreams,nUE,'uint16');
UE_used_CQI            = zeros(maxstreams,nUE,'uint8');

UE_ACK_codeblocks      = zeros(maxstreams,nUE,'uint16');
UE_C                   = zeros(maxstreams,nUE,'uint8');
UE_avg_CB_size         = zeros(maxstreams,nUE);

% Loop over all UEs and streams
for uu = 1:nUE
    for stream_i = 1:BS_output.UE_signaling(uu).MCS_and_scheduling.nCodewords
        
        % Update ACK (UE traces) and biterrors (cell traces)
        UE_ACK(stream_i,uu) = UE_output(uu).ACK(stream_i);
        UE_ACK_codeblocks(stream_i,uu) = UE_output(uu).ACK_codeblocks(stream_i);
        UE_C(stream_i,uu) = UE_output(uu).C(stream_i);
        UE_avg_CB_size(stream_i,uu) = mean(BS_output.UE_signaling(uu).TB_segmentation(stream_i).CB_sizes);
        
        UE_rv_idx(stream_i,uu) = UE_output(uu).rv_idx(stream_i);
        UE_RBs_assigned(uu) = BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs;
        
        % Update biterrors (UE traces)
        %                     if(length(UE_output(uu).rx_data_bits{stream_i}) > length(BS_output.genie(uu).data_bits{stream_i}))
        %                         help  = UE_output(uu).rx_data_bits{stream_i};
        %                         UE_output(uu).rx_data_bits{stream_i} = help(1:length(BS_output.genie(uu).data_bits{stream_i}));
        %                     elseif(length(UE_output(uu).rx_data_bits{stream_i}) < length(BS_output.genie(uu).data_bits{stream_i}))
        %                         help = BS_output.genie(uu).data_bits{stream_i};
        %                         BS_output.genie(uu).data_bits{stream_i} = help(1:length(UE_output(uu).rx_data_bits{stream_i}));
        %                     end
        
        if UE_output(uu).UE_scheduled
            
            UE_biterrors_coded(stream_i,uu)   = sum(abs(UE_output(uu).rx_data_bits{stream_i}  - BS_output.genie(uu).data_bits{stream_i}));
            UE_biterrors_uncoded(stream_i,uu) = sum(abs(UE_output(uu).rx_coded_bits{stream_i} - BS_output.genie(uu).sent_bits{stream_i}));
            
            % Update blocksize (UE traces)
            UE_blocksize_coded(stream_i,uu)   = length(UE_output(uu).rx_data_bits{stream_i});
            UE_blocksize_uncoded(stream_i,uu) = length(UE_output(uu).rx_coded_bits{stream_i});
            
            % Update FER and throughput (UE traces)
            % Coded
            if UE_output(uu).ACK(stream_i)
                UE_throughput_coded(stream_i,uu) = length(UE_output(uu).rx_data_bits{stream_i});
                UE_throughput_useful(stream_i,uu) = BS_output.UE_signaling(uu).MCS_and_scheduling.N_used_bits(stream_i);
                UE_FER_coded(stream_i,uu) = 0;
            else
                UE_FER_coded(stream_i,uu) = 1;
                UE_throughput_coded(stream_i,uu) = 0;
                UE_throughput_useful(stream_i,uu) = 0;
            end
            
            % Uncoded
            if (UE_biterrors_uncoded(stream_i)==0)
                UE_throughput_uncoded(stream_i,uu) = length(UE_output(uu).rx_coded_bits{stream_i});
                UE_FER_uncoded(stream_i,uu) = 0;
            else
                UE_FER_uncoded(stream_i,uu) = 1;
                UE_throughput_uncoded(stream_i,uu) = 0;
            end
            
        end
        
        % Update what codewords were used (valid positions in the traces, 0 if no RBs were allocated)
        if(BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs)
            UE_used_codewords(stream_i,uu) = 1;
            UE_used_CQI(stream_i,uu) = BS_output.UE_signaling(uu).MCS_and_scheduling.cqi(stream_i);
        else
            UE_used_codewords(stream_i,uu) = 0;
        end
        
    end
    
    % Update cell coded and uncoded bit errors
    cell_biterrors_coded(:) = cell_biterrors_coded(:) + UE_biterrors_coded(:,uu);
    cell_biterrors_uncoded(:) = cell_biterrors_uncoded(:) + UE_biterrors_uncoded(:,uu);
    
    % Update blocksize (cell traces)
    cell_blocksize_coded(:)   = cell_blocksize_coded(:)  + UE_blocksize_coded(:,uu);
    cell_blocksize_uncoded(:) = cell_blocksize_uncoded(:) + UE_blocksize_uncoded(:,uu);
    
    % Update FER and throughput (cell traces)
    cell_throughput_coded(:)  = cell_throughput_coded(:) + UE_throughput_coded(:,uu);
    cell_throughput_uncoded(:)  = cell_throughput_uncoded(:) + UE_throughput_uncoded(:,uu);
    cell_throughput_useful(:)  = cell_throughput_useful(:) + UE_throughput_useful(:,uu);
    cell_FER_coded(:) = uint16(cell_FER_coded(:))+ uint16(UE_FER_coded(:,uu));
    cell_FER_uncoded(:) = uint16(cell_FER_uncoded(:))+ uint16(UE_FER_uncoded(:,uu));
    
    % The number of codewords is the maximum (bitwise OR) of the used codewords matrix
    cell_used_codewords(:) = uint16(cell_used_codewords(:)) + uint16(UE_used_codewords(:,uu));
    
    % If 0 RBs were allocated => no transmission  -> no channel error calculation for this user
    if(BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs)
        channel_error_help = zeros(size(cell_channel_error));
        channel_error_help = UE_output(uu).channel_estimation_error;
        cell_channel_error = cell_channel_error + channel_error_help;
    end
end
cell_channel_error = cell_channel_error/nUE;

end
        