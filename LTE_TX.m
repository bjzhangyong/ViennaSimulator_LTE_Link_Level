function LTE_TX(LTE_params, BS, UE, AtPort, subframe_i, UE_output,BS_output)
% LTE transmitter.
% [BS_output] = LTE_TX(BS, UE, AtPort, subframe_i)
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   BS                  ... [1 x nBS]struct  - Base stations parameters from LTE_load_parameters.m
%           UE                  ... [1 x nUE]struct  - User equipments capabilities from LTE_load_parameters.m
%           AtPort              ... [1 x 1]double  - antenna port number
%           subframe_i          ... [1 x 1]double  - number of the subframe transmitted
%           previous_ACK        ... [1 x 1]logical - whether the last ACK
%                                                    was correct or not
%           BS_output           ... [1 x nBS]struct - Base Station output
%  output:  UE                  ... [1 x nUE]struct  - User equipments capabilities from LTE_load_parameters.m
%
% date of creation: 2008/08/11
% last changes: 2008/09/02  Colom Ikuno HARQ is implemented
%               2008/09/08  Bosanska    Rewritten for the multi-user scenario
%               2008/09/10  Bosanska    previous_ACK is now stored in
%                                       UE_output [1 x nUE] struct for UE output (ACK for multiple users)
%               2008/09/15  Bosanska    added updated struct BS as an output
%                                       changed structure of BS and UE (acc. to LTE_load_parameters.m)
%                                       changed structure of BS_output
%               2008/10/21  Bosanska    added updated struct BS as an output
%                                       changed structure of BS and UE (acc. to LTE_load_parameters.m)
%                                       changed structure of BS_output

if LTE_params.UE_config.mode == 6 % Interference Alignment
    V = BS_output.UE_signaling.MCS_and_scheduling.V;
    U = BS_output.UE_signaling.MCS_and_scheduling.U;
end

%% Calculate the correct number of subframe from 1 - 10
subframe_corr = mod(subframe_i,10);
if(subframe_corr == 0)
    subframe_corr = 10;
end

%% Read pregenerated reference signals

RefSym          = LTE_params.Reference_Signal(1,subframe_corr).RefSym;
RefMapping      = LTE_params.Reference_Signal(1,subframe_corr).RefMapping;
total_no_refsym = LTE_params.Reference_Signal(1,subframe_corr).total_no_refsym;     % total number of reference symbols per slot per resource block
NoData_indices  = LTE_params.Reference_Signal(1,subframe_corr).NoData_indices;      % indices of resource elements where no data is to be written (for the whole subframe)
% [RefSym,RefMapping,total_no_refsym,NoData_indices]=LTE_Common_gen_Reference_Signal(LTE_params.Nrb, LTE_params.Nsub, BS.nTX, BS.NIDcell, subframe_corr);

%% Generation of other channels (PBCH, PDSCH)
CHmapping = zeros(size(NoData_indices));
CHusedElements = zeros(LTE_params.Nrb*2,1);
CHnsyms_subframe = 0;
if LTE_params.usePBCH
    CHmapping = CHmapping | LTE_params.PBCH(1,subframe_corr).Mapping;
    CHusedElements = CHusedElements + LTE_params.PBCH(1,subframe_corr).UsedElements;
%     CHnsyms_subframe = CHnsyms_subframe + LTE_params.PBCH(1,subframe_corr).N_Elements;
end
if LTE_params.usePDCCH
    CHmapping = CHmapping | LTE_params.PDCCH(1,subframe_corr).Mapping;
    CHusedElements = CHusedElements + LTE_params.PDCCH(1,subframe_corr).UsedElements;
%     CHnsyms_subframe = CHnsyms_subframe + LTE_params.PDCCH(1,subframe_corr).N_Elements;
end

%% Generation of synchronization signals and Initialization of signals
PrimMapping = zeros(size(NoData_indices));
SecMapping = zeros(size(NoData_indices));
PrimSync = [];
SecSync = [];
if(subframe_corr == 1 || subframe_corr == 6)
    % Read pregenerated sync signals
    PrimSync         = LTE_params.Sync_Signal(1,subframe_corr).PrimSync;
    PrimMapping      = LTE_params.Sync_Signal(1,subframe_corr).PrimMapping;
    SecSync          = LTE_params.Sync_Signal(1,subframe_corr).SecSync;
    SecMapping       = LTE_params.Sync_Signal(1,subframe_corr).SecMapping;
    SyncUsedElements = LTE_params.Sync_Signal(1,subframe_corr).SyncUsedElements;
    ResMapping = LTE_params.Sync_Signal(1,subframe_corr).ResMapping;    % reserved REs, nothing is transmitted here!
    CHusedElements = CHusedElements + LTE_params.Sync_Signal(1,subframe_corr).ResUsedElements;
%     CHnsyms_subframe = CHnsyms_subframe + LTE_params.Sync_Signal(1,subframe_corr).ResNElements;
    CHmapping = CHmapping | ResMapping;
    % [PrimSync, PrimMapping, SecSync, SecMapping, SyncUsedElements] = LTE_Common_gen_Synchronization_Signal(BS.NIDcell, LTE_params.Nrb, LTE_params.Nsc, LTE_params.Nsub, subframe_corr);
else
    for ii = 1:BS.nAtPort
        %tx_symbols{ii} = zeros(LTE_params.Nrb,2,LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))));    
                tx_symbols{ii} = zeros((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
    end
    SyncUsedElements = zeros(LTE_params.Nrb,2);
end
y_tx_assembled = zeros(LTE_params.Ntot,2*LTE_params.Ns,BS.nTX);
if LTE_params.usePBCH && subframe_corr == 6
    CHusedElements(:) = LTE_params.Ns*LTE_params.Nsc;
    CHusedElements = CHusedElements - SyncUsedElements(:);
    CHusedElements = CHusedElements -total_no_refsym;

%     CHmapping = ones(size(CHmapping))- (PrimMapping | SecMapping);
end
% image((RefMapping | PrimMapping | SecMapping |CHmapping)*100)
% grid on
% gridcolor(gca,'w','w','w')
% set(gca,'Fontsize',24,'Linewidth',1.5,'XTick',[0.5:7:14.5],'YTick',[0.5:12:72.5],'YTicklabel',[1,13,25,37,49,61,73]-1,'XTicklabel',[1,8,15]-1)
% xlabel('Sample number')
% ylabel('Subcarrier number')

%% Scheduling
% The scheduler now calculates the number of allocated bits. This is why it needs the number of reference symbols.
UE_MCS_and_scheduling_info = BS.scheduler.scheduler_users(subframe_corr,total_no_refsym,SyncUsedElements,UE_output,BS.UE_specific,BS_output.cell_genie,CHusedElements);
% UE_MCS_and_scheduling_info = BS.scheduler.scheduler_users(subframe_corr,total_no_refsym,SyncUsedElements,UE_output,BS.UE_specific,BS_output.cell_genie);
%% Transmission for each user
for uu = 1:LTE_params.nUE
    % To avoid problems with old information making it to the next loop iteration
    BS_output.UE_signaling(uu).clear;
    
    % Add the scheduler information to the downlink channel
    BS_output.UE_signaling(uu).MCS_and_scheduling = UE_MCS_and_scheduling_info(uu);
    
    if LTE_params.UE_config.mode == 6 % Interference Alignment
        BS_output.UE_signaling(uu).MCS_and_scheduling.V = V;
        BS_output.UE_signaling(uu).MCS_and_scheduling.U = U;
    end
    
    % Get necessary data and convert in in the correct format
    nLayers      = UE_MCS_and_scheduling_info(uu).nLayers;
    nCodewords   = UE_MCS_and_scheduling_info(uu).nCodewords;
    tx_mode      = UE_MCS_and_scheduling_info(uu).tx_mode;
    UE_mapping   = UE_MCS_and_scheduling_info(uu).UE_mapping;
    assigned_RBs = UE_MCS_and_scheduling_info(uu).assigned_RBs;

    if(BS_output.UE_signaling(uu).MCS_and_scheduling.assigned_RBs)

        % Get what HARQ process was assigned to this TX
        UE_output(uu).HARQ_process = BS.UE_specific(uu).current_HARQ_process.id;

        layer_x = []; % clear layer_x for every user
        for cw = 1:nCodewords

            % Get the number of data and coded bits from the scheduler
            N_coded_bits = UE_MCS_and_scheduling_info(uu).N_coded_bits(cw);
            N_data_bits  = UE_MCS_and_scheduling_info(uu).N_data_bits(cw);
            tx_rv_idx = BS.UE_specific(uu).current_HARQ_process(cw).rv_idx;
            % Generate data bits
            if ~strcmp(LTE_params.Simulation_type,'TB')
                if tx_rv_idx==0
                    if LTE_params.simulate_with_all_zero_sequences
                        tx_data_bits = false(1,N_data_bits);
                    else
                        tx_data_bits = logical(randi(LTE_params.data_RandStream,[0,1],[1,N_data_bits]));
                    end
                    BS.UE_specific(uu).current_HARQ_process(cw).HARQ_tx_buffer = tx_data_bits;
                else
                    % Retransmission, use previously stored bits in the current HARQ process
                    tx_data_bits = BS.UE_specific(uu).current_HARQ_process(cw).HARQ_tx_buffer;
                end
            else
                if LTE_params.simulate_with_all_zero_sequences
                    tx_data_bits = false(1,N_data_bits);
                else
                    tx_data_bits = logical(randi(LTE_params.data_RandStream,[0,1],[1,N_data_bits]));
                end
                BS.UE_specific(uu).current_HARQ_process(cw).HARQ_tx_buffer = tx_data_bits;
            end
            %% Coding of the bits, as of TS 36.212
            BS_output.UE_signaling(uu).turbo_rate_matcher(cw).G = N_coded_bits; % How many bits the we are allowed to transmit. Decided by the scheduler(uu)
            % According to TS 36.212, subclause 5.1.4.1.2.
            % NL is equal to 1 for blocks mapped onto one transmission layer and
            % is equal to 2 for blocks mapped onto two or four transmission layers
%             codeword_layers = nLayers;    % layers used for the actual codeword
            switch nLayers
                case 1
                    BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 1;
                case 2
                    BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 2;
                case 3
                    if nCodewords == 1
                        BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 3;   % NOTE: this might be wrong: 36.212 does not define what happens if we have 3 layers per codeword!
                    else
                        if (cw == 1)
                            BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 1; 
                        else
                            BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 2;
                        end
                    end
                case 4
                    BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 2;
%                 case 5
%                     if cw == 1
%                         BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 2;   
%                     else
%                         BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 3;   % NOTE: this might be wrong: 36.212 does not define what happens if we have 3 layers per codeword!
%                     end
%                 case 6
%                         BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 3;   % NOTE: this might be wrong: 36.212 does not define what happens if we have 3 layers per codeword!
%                 case 7
%                     if cw == 1
%                         BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 3;   % NOTE: this might be wrong: 36.212 does not define what happens if we have 3 layers per codeword!
%                     else
%                         BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 2;
%                     end
%                 case 8 
%                     BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = 2;
            end

            %BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l = min(BS_output.UE_signaling(uu).MCS_and_scheduling.nLayers,BS.nTX); % This should be a new variable "layer number" or something like that % 1=map to 1 transmission layer, 2 for 2 or 4 transmission layers
%             BS_output.UE_signaling(uu).turbo_rate_matcher(cw).N_l    = min(codeword_layers,BS.nTX);
            BS_output.UE_signaling(uu).turbo_rate_matcher(cw).rv_idx = tx_rv_idx; % redundancy version index

            tx_coded_bits = LTE_tx_DLSCH_encode(LTE_params,tx_data_bits,BS_output.UE_signaling(uu),BS_output.genie(uu),UE(uu),cw);

            %% Scrambling of the bits, as of TS 36.211 V8.2.0 (2008-03) Section 6.3.1
            tx_scrambled_bits = LTE_common_scrambling(tx_coded_bits,BS.NIDcell,uu,subframe_corr,cw,'scramble');

            %% Symbol mapping
            % NOTE: why is this called 'nibble'?
            switch BS_output.UE_signaling(uu).MCS_and_scheduling.CQI_params(cw).modulation_order
                case 1
                    nibble = tx_scrambled_bits;
                case 2
                    nibble = 2 * tx_scrambled_bits(1:2:end) + tx_scrambled_bits(2:2:end);
                case 4
                    nibble = 8 * tx_scrambled_bits(1:4:end) + 4 * tx_scrambled_bits(2:4:end) + 2 * tx_scrambled_bits(3:4:end) + tx_scrambled_bits(4:4:end);
                case 6
                    nibble = 32* tx_scrambled_bits(1:6:end) + 16* tx_scrambled_bits(2:6:end) + 8 * tx_scrambled_bits(3:6:end) + 4 * tx_scrambled_bits(4:6:end) + 2 * tx_scrambled_bits(5:6:end) + tx_scrambled_bits(6:6:end);
                otherwise
                    error('Mapping not implemented  for this constellation');
            end
            if mod(length(nibble),2) ~= 0
                nibble = nibble(1:end-1);
            end
            tx_user_symbols = LTE_params.SymbolAlphabet{BS_output.UE_signaling(uu).MCS_and_scheduling.CQI_params(cw).modulation_order}(nibble+1).';

            %% Layer mapping
            layer_x = LTE_common_layer_mapping(tx_mode,tx_user_symbols,nLayers,layer_x,nCodewords,cw,BS.nAtPort);
        end
  
        %% Determination of Codebook index
        if (tx_mode == 4) % closed loop spatial multiplexing
            if (BS_output.UE_signaling(uu).MCS_and_scheduling.CDD == 2) % CDD = 2 corresponds large delay CDD
                error('Only small or zero delay CDD is supported for closed loop spatial multiplexing');
            end
            BS_output.UE_signaling(uu).MCS_and_scheduling.codebook_index = BS_output.UE_signaling(uu).MCS_and_scheduling.PMI;
        else
            BS_output.UE_signaling(uu).MCS_and_scheduling.codebook_index = 0;
        end

        if (tx_mode == 3)    % open loop Spatial multiplexing
%             if (UE_output(uu).RI == 1) % use Tx-Diversity for RI = 1 according to 3GPP TS 36.213-820, Section 7.1.3 page 13 % NOTE: maybe this should be checked in the scheduler
%                 BS_output.UE_signaling(uu).MCS_and_scheduling.tx_mode = 2; 
%             else BS_output.UE_signaling(uu).MCS_and_scheduling.CDD = 2; % use large delay CDD
%             end
            BS_output.UE_signaling(uu).MCS_and_scheduling.CDD = 2; % use large delay CDD
            if (BS.nAtPort == 2) % choice of precoding matrix according to 3GPP TS 36.213-820, Section 7.1.3 page 13
                BS_output.UE_signaling(uu).MCS_and_scheduling.codebook_index = 0;
            end
            if (BS.nAtPort == 4)
                BS_output.UE_signaling(uu).MCS_and_scheduling.codebook_index = [12,13,14,15];
            end
            if (BS.nAtPort == 8)
                error('still missing');
            end
        end

        %% Pre-Coding
        
        RB_indices = find(UE_mapping == 1);
        BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices = [];
        BS_output.UE_signaling(uu).MCS_and_scheduling.slot_indices = [];
%         assigned_RBs
        for ww = 1:assigned_RBs
            zero_temp = zeros(size(UE_mapping));
            zero_temp(RB_indices(ww)) = 1;
            if(subframe_corr == 1 || subframe_corr == 6)
                freq_tmp = mod(find(kron((zero_temp.*UE_mapping),ones(LTE_params.Nsc,LTE_params.Ns)).*~(NoData_indices|PrimMapping|SecMapping|CHmapping) ~= 0) ,LTE_params.Ntot);
                BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices = [BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices;freq_tmp];
                BS_output.UE_signaling(uu).MCS_and_scheduling.slot_indices = [BS_output.UE_signaling(uu).MCS_and_scheduling.slot_indices;(floor((RB_indices(ww)-1)/LTE_params.Nrb)+1)*ones(size(freq_tmp))];
            else
                freq_tmp = mod(find(kron((zero_temp.*UE_mapping),ones(LTE_params.Nsc,LTE_params.Ns)).*~(NoData_indices|CHmapping) ~= 0) ,LTE_params.Ntot);
                BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices = [BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices;freq_tmp];
                BS_output.UE_signaling(uu).MCS_and_scheduling.slot_indices = [BS_output.UE_signaling(uu).MCS_and_scheduling.slot_indices;(floor((RB_indices(ww)-1)/LTE_params.Nrb)+1)*ones(size(freq_tmp))];
            end
        end
        BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices(~BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices) = LTE_params.Ntot; % where the mod operation delivers a zero there should be Ntot
        % This is the actual precoding part
        [precode_y,BS_output.UE_signaling(uu).MCS_and_scheduling.PRE] = LTE_precoding(...
            BS_output.UE_signaling(uu).MCS_and_scheduling.tx_mode,...
            layer_x,...
            BS.nAtPort,...
            BS_output.UE_signaling(uu).MCS_and_scheduling.codebook_index,...
            BS_output.UE_signaling(uu).MCS_and_scheduling.nLayers,...
            BS_output.UE_signaling(uu).MCS_and_scheduling.freq_indices,...
            LTE_params.Ntot,...
            BS_output.UE_signaling(uu).MCS_and_scheduling.CDD,...
            LTE_params,...
            BS_output.UE_signaling(uu).MCS_and_scheduling.slot_indices,...
            BS_output.UE_signaling(uu).MCS_and_scheduling.V,...
            BS_output.UE_signaling(uu).MCS_and_scheduling.U);

        BS_output.genie(uu).layer_x = layer_x;
        %% OFDM symbol assembly
        % Based on 3GPP TS 36.211 V.8.2.0 section 6.3.5
        % The mapping to resource elements (k,l) on antenna port p not reserved for other purposes in increasing order
        % of first the index k, and then the index l, starting with the first slot in a subframe.

        for nn = 1:BS.nAtPort
            rb_numbers = find(UE_mapping == 1);
            if(subframe_corr == 1 || subframe_corr == 6 || LTE_params.usePDCCH)
                rb_tx_symbols = LTE_params.Nsc*LTE_params.Ns - total_no_refsym - SyncUsedElements(UE_mapping) - CHusedElements(UE_mapping);
                index = 0;
                for ii = 1:length(rb_numbers)
                    if(rb_numbers(ii) > LTE_params.Nrb)
                        y_tx_assembled_temp = zeros(LTE_params.Nsc,LTE_params.Ns); 
                        y_tx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)+CHmapping((rb_numbers(ii)-1-LTE_params.Nrb)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,8:8+LTE_params.Ns-1))) = precode_y(nn,index+1:index+rb_tx_symbols(ii)).';
                        y_tx_assembled((rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end,nn) = y_tx_assembled_temp;
                    else
                        y_tx_assembled_temp = zeros(LTE_params.Nsc,LTE_params.Ns);
                        y_tx_assembled_temp(~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)...
                            + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)+...
                            CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns))) = precode_y(nn,index+1:index+rb_tx_symbols(ii)).';
                        y_tx_assembled((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns,nn) = y_tx_assembled_temp;
                    end
                    index = index + rb_tx_symbols(ii);
                end
            else
                % NOTE: this may give problems. What UE_mapping? There are 2,
                % one for each codeword!!!! - i'm not sure about this...
                tx_symbols{nn}(repmat(reshape(UE_mapping,[1 size(UE_mapping)]),[LTE_params.Nsc*LTE_params.Ns - ...
                    sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))),1,1])) = precode_y(nn,:);
                %tx_symbols(:,:,:,nn)=tx_symbols_temp{nn};
            end
        end
    else
        %if specific user is not scheduled, then the genie data from previous frame have to be deleted
        for cw = 1:length(BS_output.genie(uu).data_bits)
            BS_output.genie(uu).data_bits{cw} = [];
            BS_output.genie(uu).sent_bits{cw} = [];
            BS_output.genie(uu).bits_to_turboencode{cw} = [];
            BS_output.cell_genie.layer_x{cw} = [];
        end
    end
end

for nn = 1:BS.nAtPort
    if(subframe_corr == 1 || subframe_corr == 6 || LTE_params.usePDCCH)
        y_tx_assembled_temp = y_tx_assembled(:,:,nn);
        y_tx_assembled_temp(logical(PrimMapping)) = PrimSync;
        y_tx_assembled_temp(logical(SecMapping)) = SecSync;
        y_tx_assembled(:,:,nn) = y_tx_assembled_temp;
    else
        y_tx_assembled_temp = zeros(LTE_params.Nsc,LTE_params.Ns);
        for ii = 1:LTE_params.Nrb
            y_tx_assembled_temp(~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)) = tx_symbols{nn}(:,ii,1);
            y_tx_assembled((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,1:LTE_params.Ns,nn) = y_tx_assembled_temp;
            y_tx_assembled_temp(~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)) = tx_symbols{nn}(:,ii,2);
            y_tx_assembled((ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,LTE_params.Ns+1:end,nn) = y_tx_assembled_temp;
        end
    end
%     y_tx_assembled(RefMapping) = RefSym(RefSym~=0);

    %% Pilot symbol boosting
    if ~LTE_params.Pilots_power_offset==0
        y_tx_assembled(:,:,nn) = LTE_params.data_factor*y_tx_assembled(:,:,nn);
        y_tx_assembled(RefMapping) = RefSym(RefSym~=0);
    else
        y_tx_assembled(RefMapping) = RefSym(RefSym~=0); 
    end
    
    %% predistortion of TB
        if strcmp(LTE_params.Simulation_type,'TB')
            %load('D:\Matlab\LTE\trunk - Copy (2)\predistortion\H_est_10_10.mat');
            H_est = LTE_params.Testbed_settings.predistortion_matrix(:,LTE_params.Testbed_settings.ant_permutations,LTE_params.Testbed_settings.Testbed_ID);
            predistortion = 1./H_est(:,nn);
            predistortion = predistortion./sqrt(mean(abs(predistortion).^2));
            y_tx_assembled(:,:,nn) = repmat(predistortion,[1 LTE_params.Nsub]).*y_tx_assembled(:,:,nn);
            %y_tx_assembled(RefMapping) = RefSym(RefSym~=0);
        end
    
    %% Zero padding, IFFT, add CP
    % insert zero DC carrier
    y_tx_assembled_padded = [y_tx_assembled(1:LTE_params.Ntot/2,:,nn); zeros(1,LTE_params.Nsub); y_tx_assembled(LTE_params.Ntot/2+1:end,:,nn)];
    % y_tx_assembled_padded2 = padarray(y_tx_assembled_padded,LTE_params.Nfft-LTE_params.Ntot-1,'post');
    % Avoid using the padarray function just for this, as it is from the signal processing toolbox
    y_tx_assembled_padded2 = [y_tx_assembled_padded; zeros(LTE_params.Nfft-LTE_params.Ntot-1,size(y_tx_assembled_padded,2))];
    y_tx_assembled_shifted = circshift(y_tx_assembled_padded2,-LTE_params.Ntot/2);

    y_tx_assembled_ifft = sqrt(LTE_params.Nfft)*sqrt(LTE_params.Nfft/LTE_params.Ntot)*ifft(y_tx_assembled_shifted);
    if(length(LTE_params.Tg)==2)
        BS_output.y_tx(:,nn) = [        y_tx_assembled_ifft(LTE_params.Index_TxCyclicPrefix{1},1);...
            reshape(y_tx_assembled_ifft(LTE_params.Index_TxCyclicPrefix{2},2:LTE_params.Ns),[],1);...
            y_tx_assembled_ifft(LTE_params.Index_TxCyclicPrefix{1},LTE_params.Ns+1);...
            reshape(y_tx_assembled_ifft(LTE_params.Index_TxCyclicPrefix{2},LTE_params.Ns+2:end),[],1) ];
    else
        BS_output.y_tx(:,nn) = reshape(y_tx_assembled_ifft(LTE_params.Index_TxCyclicPrefix,:),[],1);
    end
end

%% mapping of antenna port to antenna
%BS_output.y_tx = zeros(size(BS_output.y_tx,2),BS.nTX);
BS_output.y_tx = LTE_map2antenna(BS_output.y_tx,AtPort+1); % now this maps antenna port {0,1,2,3} to antennas {1,2,3,4}

%% Fill in extra genie information
BS_output.cell_genie.y_tx_assembled = y_tx_assembled;

