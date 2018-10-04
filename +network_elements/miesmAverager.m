classdef miesmAverager < network_elements.sinrAverager
    % Implements an MIESM averager. Given a set of predefined \beta values
    % and the Modulation and Coding Scheme (MCS) used, returns an effective
    % SINR value.
    % (c) Stefan Schwarz 2010


   properties
       % Place where we will store the beta value for each possible MCS
       beta_tables
       % MCSs corresponding to the previous beta values
       MCS_values
       % Indicate which MCSs are valid
       isvalid
       % BICM data
       MI_data
   end

   methods
       function obj = miesmAverager(beta_values,MCS_values)
           % Input checking
           if length(beta_values) ~= length(MCS_values)
               error('The vector containg the beta values and the vector specifying the corresponging MCSs are of different length');
           end
           if min(MCS_values)<0
               error('The minimum MCS cannot be lower than 0');
           end
           
           % Translate the MCSs so we can put them in a lookup table.
           for i_=1:length(beta_values)
               MCS_value_idx = MCS_values(i_) + 1;
               obj.MCS_values(MCS_value_idx)  = MCS_values(i_);
               obj.beta_tables(MCS_value_idx) = beta_values(i_);
               obj.isvalid(MCS_value_idx)     = true;
           end
           MI_data_loader = utils.MIdata_loader;
           obj.MI_data = MI_data_loader.MIdata_load;
       end
       
       % Average the given SINR vector. varargin contains the following:
       %   - varargin{1} = MCS -> values in the range 0:15
       %   - varargin{2} = symbol alphabet order corresponding to the MCS values
       % INPUT VALUES IN LINEAR
       % OUTPUT VALUE IN dB
       % Allows for calculations with multiple beta values
       function effective_SINR = average(obj,SINR_vector,varargin)
           % Ensure matrix dimensions
           SINR_vector = SINR_vector(:);
           MCS_idx = varargin{1};
           alphabets = varargin{2}/2;
           effective_SINR = zeros(length(MCS_idx),1);
           if min(MCS_idx)<0
               error('CQI cannot be lower than 0');
           elseif max(MCS_idx)>length(obj.isvalid)
               error('CQI cannot be higher than %d',length(obj.isvalid));
           elseif sum(~obj.isvalid(MCS_idx+1))~=0
               error('SINR averaging not defined for all of the CQIs');
           end
           if length(SINR_vector) > 1
               for cqi_i = MCS_idx
                   %                if cqi_i
                   SINR_vector_temp = SINR_vector/obj.beta_tables(cqi_i+1);
                   %                     SINR_vector_temp(SINR_vector_temp > obj.MI_data(alphabets(cqi_i+1)).max_SNR) = obj.MI_data(alphabets(cqi_i+1)).max_SNR;
                   %                     SINR_vector_temp(SINR_vector_temp < obj.MI_data(alphabets(cqi_i+1)).min_SNR) = obj.MI_data(alphabets(cqi_i+1)).min_SNR;
                   for SNR_i = 1:length(SINR_vector_temp)
                       [~,ind(SNR_i)] = min(abs(SINR_vector_temp(SNR_i)-obj.MI_data(alphabets(cqi_i+1)).SNR));
                   end
                   I_avg = mean(obj.MI_data(alphabets(cqi_i+1)).BICM(ind));
                   [null_var,ind2] = min(abs(I_avg-obj.MI_data(alphabets(cqi_i+1)).BICM));
                   effective_SINR(cqi_i+1) = obj.beta_tables(cqi_i+1)*obj.MI_data(alphabets(cqi_i+1)).SNR(ind2);
                   %                 I_avg = mean(interp1(obj.MI_data(alphabets(cqi_i)).SNR,obj.MI_data(alphabets(cqi_i)).BICM,SINR_vector_temp,'linear','extrap'));
                   %                 effective_SINR(cqi_i+1) = obj.beta_tables(cqi_i)*interp1(obj.MI_data(alphabets(cqi_i)).BICM(1:obj.MI_data(alphabets(cqi_i)).sat_SNR),obj.MI_data(alphabets(cqi_i)).SNR(1:obj.MI_data(alphabets(cqi_i)).sat_SNR),I_avg,'linear','extrap');
                   %                else
                   %                    effective_SINR(cqi_i+1) = 10^-10;
                   %                end
               end
           else
               effective_SINR = SINR_vector;
           end
           effective_SINR = 10*log10(effective_SINR);
       end
   end
end
