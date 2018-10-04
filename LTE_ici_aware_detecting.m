function [LLR_SD M rx_layer_x_equalized interference_power signal_power] = LTE_ici_aware_detecting(LTE_params,UE_signaling,rx_user_symbols,subframe_corr,ChanMod,rb_numbers,y_rx_assembled,layer_x,H_fft_matrix_est,uu,RefSym,RefMapping,NoData_indices,varargin)
%function [LLR_SD M rx_layer_x_equalized interference_power signal_power] = LTE_ici_aware_detecting(LTE_params,UE_signaling,rx_user_symbols,subframe_corr,ChanMod,rb_numbers,y_rx_assembled,layer_x,H_fft_matrix_est,uu,rb_rx_symbols,RefSym,RefMapping,NoData_indices,PrimMapping,SecMapping,PrimSync,SecSync,CHmapping)
% LTE channel estimator - to filter the output of the transmitter.
% [chan_output] = LTE_channel_model(BS_output, SNR)
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) by INTHFT
% www.nt.tuwien.ac.at
%
% input :   
% output:   
%
% date of creation: 2011/01/27
% last changes: 2011/01/27  Simko

            nr_inter = LTE_params.subcarrier_intervals;   %number of subcarriers from which the ICI is considered
            nLayers = UE_signaling.MCS_and_scheduling.nLayers;
            rx_layer_x_equalized = nan(length(rx_user_symbols),nLayers);
            rx_layer_x_equalized_temp = nan(nLayers,LTE_params.Ntot,LTE_params.Nsub);
            
            optargin = size(varargin,2);
            if optargin==6
                PrimMapping = varargin{1};
                SecMapping = varargin{2};
                PrimSync = varargin{3};
                SecSync = varargin{4};
                CHmapping = varargin{5};
                rb_rx_symbols = varargin{6};
            end
            
            
            M = [UE_signaling.MCS_and_scheduling.CQI_params.modulation_order];
            switch nLayers % when layer number is unequal to codeword number we need to do something
                case 2
                    if(UE_signaling.MCS_and_scheduling.nCodewords == 1) % 1 codewords, 2 layers
                        M = [M,M];
                    end
                case 3  % 2 codewords, 3 layers
                    M = [M(1),M(2),M(2)];
                case 4  % 2 codewords, 4 layers
                    M = [M(1),M(1),M(2),M(2)];
            end
            bittable = false(sum(M(1:nLayers)),2^max(M));
            symbol_alphabet = zeros(nLayers,2^max(M));
            for i = 1:nLayers
                bittable(sum(M(1:i-1))+(1:M(i)),1:2^M(i))=LTE_params.bittable{M(i)}; % Bitmapping table
                symbol_alphabet(i,1:2^M(i))=LTE_params.SymbolAlphabet{M(i)}.'; % Symbol alphabet
            end
            LLR_SD = zeros(sum(M),length(rx_user_symbols));   % Log likelihood Ratios of the Spere decoder
            LLR_SD_temp = zeros(sum(M),LTE_params.Ntot,LTE_params.Nsub);
            
            if(subframe_corr == 1 || subframe_corr == 6)
                ref_symbols = zeros(LTE_params.Ntot,LTE_params.Nsub,ChanMod.nTX);
                ref_symbols(RefMapping) =  RefSym(RefSym~=0);
                ref_symbols(repmat(PrimMapping,[1,1,ChanMod.nTX])) = repmat(PrimSync,1,ChanMod.nTX);
                ref_symbols(repmat(SecMapping,[1,1,ChanMod.nTX])) = repmat(SecSync,1,ChanMod.nTX);
                no_datas = PrimMapping | SecMapping | NoData_indices | CHmapping;
                
            else
                ref_symbols = zeros(LTE_params.Ntot,LTE_params.Nsub,ChanMod.nTX);
                ref_symbols(RefMapping) =  RefSym(RefSym~=0);
                no_datas = NoData_indices;
                
            end
            
            precoding_used = zeros(LTE_params.Ntot,LTE_params.Nsub);
            precoding_nr = size(UE_signaling.MCS_and_scheduling.PRE,3);
            precoding_vec = repmat(1:precoding_nr,1,ceil(length(rx_user_symbols)/precoding_nr));
            precoding_vec = precoding_vec(1:length(rx_user_symbols)).';
            
            freq_index = 1:12;
            precoding_start = 1;
            for rb_i =1:length(rb_numbers)
                precoding_temp = zeros(LTE_params.Nsc,LTE_params.Ns);
                if rb_numbers(rb_i)>LTE_params.Nrb
                    time_index = (1:LTE_params.Ns) + 7;
                else
                    time_index = (1:LTE_params.Ns);
                end
                precoding_length = sum(sum(~no_datas(freq_index,time_index)));
                precoding_temp(~no_datas(freq_index,time_index)) = precoding_vec(precoding_start:precoding_start+precoding_length-1);
                precoding_used(freq_index,time_index) = precoding_temp;
                
                
                precoding_start = precoding_start + precoding_length;
                
                freq_index = mod(freq_index +12,LTE_params.Ntot);
                freq_index(find(freq_index==0)) = LTE_params.Ntot;
            end
            
            y_rx_symbol_ref_free = zeros(size(y_rx_assembled));
            y_rx_zero_force = zeros(LTE_params.Ntot,LTE_params.Nsub,nLayers);
            
            data_nr_total = 0;
            interference_matrix = nan(length(rx_user_symbols),1);     
            data_temp = nan(nLayers,72,14);

            index123 = 1;
            for rb_i = 1:12
                data_small_temp1 = zeros(LTE_params.Nsc,LTE_params.Ns);
                data_small_temp2 = zeros(LTE_params.Nsc,LTE_params.Ns);
                if(rb_i > LTE_params.Nrb)
                    no_data_small_temp = no_datas((rb_i-7)*LTE_params.Nsc+1:(rb_i-6)*LTE_params.Nsc,LTE_params.Ns+1:end);
                    nr_of_data_symbols_temp = sum(sum(~no_data_small_temp));
                    data_small_temp1(~no_data_small_temp) = layer_x(1,index123:index123+nr_of_data_symbols_temp-1);
                    data_temp(1,(rb_i-7)*LTE_params.Nsc+1:(rb_i-6)*LTE_params.Nsc,LTE_params.Ns+1:end) = data_small_temp1;
                    
                    if nLayers==2
                        data_small_temp2(~no_data_small_temp) = layer_x(2,index123:index123+nr_of_data_symbols_temp-1);
                        data_temp(2,(rb_i-7)*LTE_params.Nsc+1:(rb_i-6)*LTE_params.Nsc,LTE_params.Ns+1:end) = data_small_temp2;
                    end
                    
                    index123 = index123 + nr_of_data_symbols_temp;
                else
                    no_data_small_temp = no_datas((rb_i-1)*LTE_params.Nsc+1:rb_i*LTE_params.Nsc,1:LTE_params.Ns);
                    nr_of_data_symbols_temp = sum(sum(~no_data_small_temp));
                    data_small_temp1(~no_data_small_temp) = layer_x(1,index123:index123+nr_of_data_symbols_temp-1);
                    data_temp(1,(rb_i-1)*LTE_params.Nsc+1:rb_i*LTE_params.Nsc,1:LTE_params.Ns) = data_small_temp1;
                    
                    if nLayers==2
                        data_small_temp2(~no_data_small_temp) = layer_x(2,index123:index123+nr_of_data_symbols_temp-1);
                        data_temp(2,(rb_i-1)*LTE_params.Nsc+1:rb_i*LTE_params.Nsc,1:LTE_params.Ns) = data_small_temp2;
                    end
                    
                    index123 = index123 + nr_of_data_symbols_temp;
                end
            end
                      
            
            for symbol_i=1:LTE_params.Nsub
                
                no_data = no_datas(:,symbol_i);
                no_data_nr = sum(no_data);
                data_nr = sum(~no_data);
                
                
                channel_matrix_symbol = squeeze(H_fft_matrix_est(:,:,symbol_i,:,:));
                channel_temp = zeros(ChanMod.nRX*(LTE_params.Ntot-no_data_nr),ChanMod.nTX*(LTE_params.Ntot-no_data_nr));
                precoding_matrix = zeros(ChanMod.nTX*(LTE_params.Ntot-no_data_nr),(LTE_params.Ntot-no_data_nr)*nLayers);
                
                for ii = 1:LTE_params.Ntot
                    iii = sum(~no_data(1:ii));
                    if no_data(ii) == 0
                        precoding_matrix((iii-1)*ChanMod.nTX+1:iii*ChanMod.nTX,(iii-1)*nLayers+1:iii*nLayers) = squeeze(UE_signaling(uu).MCS_and_scheduling.PRE(:,:,precoding_used(ii,symbol_i)));
                    end
                end
                
                %Substraction of ICI from pilot and reference symbols
                for rr=1:ChanMod.nRX
                    y_rx_symbol = squeeze(y_rx_assembled(:,symbol_i,rr));
                    ref_interrerence = zeros(LTE_params.Ntot,1);
                    for tt=1:ChanMod.nTX
                        ref_symbol = squeeze(ref_symbols(:,symbol_i,tt));
                        
                        %channel_small_temp = squeeze(channel_matrix_symbol(:,:,rr,tt));
                        channel_small_temp_mapping = diag(true(LTE_params.Ntot,1));
                        
                        for sub_iii = 1:nr_inter
                            channel_small_temp_mapping = channel_small_temp_mapping + diag(true(LTE_params.Ntot-sub_iii,1),sub_iii);
                            channel_small_temp_mapping = channel_small_temp_mapping + diag(true(LTE_params.Ntot-sub_iii,1),-sub_iii);
                        end
                        channel_small_temp_mapping = logical(channel_small_temp_mapping);
                        %                         channel_small_temp1 = diag(diag(channel_small_temp));
                        %                         for sub_iii = 1:nr_inter
                        %                             channel_small_temp1 = channel_small_temp1 + diag(diag(channel_small_temp,sub_iii),sub_iii);
                        %                             channel_small_temp1 = channel_small_temp1 + diag(diag(channel_small_temp,-sub_iii),-sub_iii);
                        %                         end
                        
                        %channel_matrix_symbol(:,:,symbol_i,rr,tt) = channel_small_temp;
                        channel_small_temp_mapping1 = false(LTE_params.Ntot,LTE_params.Ntot,ChanMod.nRX,ChanMod.nTX);
                        channel_small_temp_mapping1(:,:,rr,tt) = channel_small_temp_mapping;
                        channel_small_temp = zeros(LTE_params.Ntot);
                        channel_small_temp(channel_small_temp_mapping) = squeeze(channel_matrix_symbol(channel_small_temp_mapping1));
                        ref_interrerence = ref_interrerence + channel_small_temp*ref_symbol;
                        
                        channel_temp(rr:ChanMod.nRX:end,tt:ChanMod.nTX:end) = channel_small_temp(~no_data,~no_data);
                        %channel_temp(rr:ChanMod.nRX:end,tt:ChanMod.nTX:end) = channel_matrix_symbol(~no_data,~no_data,rr,tt);
                    end
                    y_rx_symbol_ref_free(:,symbol_i,rr) = y_rx_symbol - ref_interrerence;
                end
                
                symbol_channel = channel_temp*precoding_matrix;
                y_rx_symbol_ref_free_temp = squeeze(y_rx_symbol_ref_free(~no_data,symbol_i,:));
                symbol_i_temp = zeros(length(y_rx_symbol_ref_free_temp)*nLayers,1);
                symbol_i_temp(1:nLayers:end)=y_rx_symbol_ref_free_temp(:,1);
                
                if nLayers==2
                    symbol_i_temp(2:nLayers:end)=y_rx_symbol_ref_free_temp(:,2);
                end
                
                
                window = zeros(LTE_params.Ntot,2*nr_inter+1);
                window(:,nr_inter+1) = [1:LTE_params.Ntot]';
                for inter_i = 1:nr_inter
                    window(:,inter_i) = [1:LTE_params.Ntot]' - (nr_inter - inter_i +1);
                    window(:,nr_inter+1+inter_i) = [1:LTE_params.Ntot]' + inter_i;
                end
                
                window_data = window(~repmat(no_data,1,2*nr_inter+1));
                window_data = reshape(window_data,[],2*nr_inter+1);
                
                window_data(find(window_data<1)) = nan;
                window_data(find(window_data>LTE_params.Ntot)) = nan;
                sub = 1;
                for ref_i = 1:LTE_params.Ntot
                    if(no_data(ref_i))
                        window_data(find(window_data==ref_i)) = nan;
                    else
                        window_data(find(window_data==ref_i)) = sub;
                        sub = sub + 1;
                    end
                end
                
                y_rx_zero_force_temp = zeros(size(symbol_i_temp));
                for sub_ii=1:LTE_params.Ntot
                    sub_i = sum(~no_data(1:sub_ii));
                    if ~no_data(sub_ii)
                        subcarriers = window_data(sub_i,:);
                        sub_down = min(subcarriers);
                        sub_up = max(subcarriers);
                        sub_of_int = (subcarriers(nr_inter+1)-sub_down)*nLayers+1:(subcarriers(nr_inter+1)-sub_down)*nLayers+nLayers;
                        
                        if nLayers==2
                            symbol_channel_small = symbol_channel(sub_down*nLayers-1:sub_up*nLayers,sub_down*nLayers-1:sub_up*nLayers);
                            inv_temp = pinv(symbol_channel_small);
                            symbol_temp = inv_temp*symbol_i_temp(sub_down*nLayers-1:sub_up*nLayers);
                        elseif nLayers==1
                            symbol_channel_small = symbol_channel(sub_down*nLayers:sub_up*nLayers,sub_down*nLayers:sub_up*nLayers);
                            inv_temp = pinv(symbol_channel_small);
                            symbol_temp = inv_temp*symbol_i_temp(sub_down*nLayers:sub_up*nLayers);
                        end
                        
                        rx_layer_x_equalized_temp(:,sub_ii,symbol_i) = symbol_temp(sub_of_int);
                        
                        %ICI measurment
                        genie_position_temp = data_nr_total + sum(~no_data(1:sub_ii));
                        genie_data_small = data_temp(:,sub_ii,symbol_i);
                        interference_matrix(genie_position_temp) = mean(abs(symbol_temp(sub_of_int) - genie_data_small).^2);
                        
                        Hg = inv_temp*symbol_channel_small;
                        noise_enhancement_tmp = sum(abs(inv_temp).^2,2);
                        noise_enhancement_tmp = noise_enhancement_tmp(sub_of_int);
                        noise_enhancement = [];
                        
                        for ii = 1:length(noise_enhancement_tmp)
                            noise_enhancement = [noise_enhancement;repmat(noise_enhancement_tmp(ii),M(i),1)];
                        end
                        
                        [C,I] = min((abs(symbol_temp*ones(1,2^M(i))-repmat(symbol_alphabet(1,:),length(symbol_temp),1)).').^2);
                        symbols_ZF = I.';    % ZF Symbols (integers)
                        
                        s_alph = symbol_alphabet(1,symbols_ZF).';
                        
                        dist_ZF = sum(abs(symbol_temp-s_alph).^2,1);   % distance to the ZF solution initial value for SS Decoder
                        % Soft Sphere Decoder
                        Hg_small = Hg(sub_of_int,sub_of_int);
                        if imag(Hg_small) == 0
                            Hg_small = complex(Hg_small);
                        end
                        
                        %LLR_SD_C = LTE_rx_soft_sd2(Hg,symbol_temp,dist_ZF,int32(symbols_ZF),int32(M),symbol_alphabet.',bittable)./noise_enhancement;
                        LLR_SD_C = LTE_rx_soft_sd2(Hg_small,symbol_temp(sub_of_int),dist_ZF,int32(symbols_ZF(sub_of_int)),int32(M),symbol_alphabet.',bittable)./noise_enhancement;
                        LLR_SD_temp(:,sub_ii,symbol_i) = LLR_SD_C;
                    end
                    
                end
                
                data_nr_total = data_nr_total + data_nr;
                
            end
            
            
            interference_power = mean(interference_matrix);
            signal_power = mean(mean(abs(layer_x).^2));
            
            
            
            
            %rx_user_symbols_ZF = zeros(sum(rb_rx_symbols),UE(uu).nRX);
            %rx_symbols_ZF  = zeros((LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
            % reorder symbols
            if(subframe_corr == 1 || subframe_corr == 6)
                index = 0;
                for ii = 1:length(rb_numbers)
                    if(rb_numbers(ii) > LTE_params.Nrb)
                        LLR_temp = LLR_SD_temp(:,(rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end);
                        LLR_SD(:,index+1:index+rb_rx_symbols(ii)) = LLR_temp(:,~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                    
                        rx_eq_temp = rx_layer_x_equalized_temp(:,(rb_numbers(ii)-LTE_params.Nrb-1)*LTE_params.Nsc+1:(rb_numbers(ii)-LTE_params.Nrb)*LTE_params.Nsc,LTE_params.Ns+1:end);
                        rx_layer_x_equalized(index+1:index+rb_rx_symbols(ii),:) = rx_eq_temp(:,~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)).';
                    
                    
                    else
                        LLR_temp = LLR_SD_temp(:,(rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns);
                        LLR_SD(:,index+1:index+rb_rx_symbols(ii)) = LLR_temp(:,~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)...
                            + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) + CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)));
                         
                        rx_eq_temp = rx_layer_x_equalized_temp(:,(rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns);
                        rx_layer_x_equalized(index+1:index+rb_rx_symbols(ii),:) = rx_eq_temp(:,~(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns) + PrimMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns)...
                            + SecMapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns) + CHmapping((rb_numbers(ii)-1)*LTE_params.Nsc+1:rb_numbers(ii)*LTE_params.Nsc,1:LTE_params.Ns))).';
                        
                    end
                    index = index + rb_rx_symbols(ii);
                end
            else
                LLR_temp2 = zeros(sum(M),(LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
                rx_eq_temp2 = zeros(nLayers,(LTE_params.Nsc*LTE_params.Ns - sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns)))),LTE_params.Nrb,2);
                for ii = 1:LTE_params.Nrb
                    % NOTE: some more comments
                    LLR_temp = LLR_SD_temp(:,(ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,1:LTE_params.Ns);
                    LLR_temp2(:,:,ii,1) = LLR_temp(:,~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                    LLR_temp = LLR_SD_temp(:,(ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,LTE_params.Ns+1:end);
                    LLR_temp2(:,:,ii,2) = LLR_temp(:,~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                    
                    rx_eq_temp = rx_layer_x_equalized_temp(:,(ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,1:LTE_params.Ns);
                    rx_eq_temp2(:,:,ii,1) = rx_eq_temp(:,~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                    rx_eq_temp = rx_layer_x_equalized_temp(:,(ii-1)*LTE_params.Nsc+1:ii*LTE_params.Nsc,LTE_params.Ns+1:end);
                    rx_eq_temp2(:,:,ii,2) = rx_eq_temp(:,~NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns));
                end
                %rx_symbols_temp = rx_symbols(:,:,:);
                LLR_SD = LLR_temp2(:,repmat(reshape(UE_signaling.MCS_and_scheduling.UE_mapping,[1 size(UE_signaling.MCS_and_scheduling.UE_mapping)]),[LTE_params.Nsc*LTE_params.Ns - ...
                    sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))),1,1]));
                
                rx_layer_x_equalized_temp2 = rx_eq_temp2(:,repmat(reshape(UE_signaling.MCS_and_scheduling.UE_mapping,[1 size(UE_signaling.MCS_and_scheduling.UE_mapping)]),[LTE_params.Nsc*LTE_params.Ns - ...
                    sum(sum(NoData_indices(1:LTE_params.Nsc,1:LTE_params.Ns))),1,1]));
                rx_layer_x_equalized = rx_layer_x_equalized_temp2.';
                
            end