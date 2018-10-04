function [ChanMod_output, SINR] = LTE_channel_matrix(LTE_params,ChanMod, SNR, UE_output, UE, uu, signal_length, subframe_i,varargin)
% LTE channel matrix - to generate channel matrix according to channel model.
% [ChanMod_output] = LTE_channel_matrix(ChanMod)
% Author: Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input :   ChanMod             ... [1x1]struct - channel model, including PDP, number of RX/TX antennas, correlation ...
% output:   ChanMod_output      ... [1x1]struct - channel output:
%                                   [nTR x nTX x length(PDP)]double channel matrix H
%
% date of creation: 2008/09/18
% last changes: 2008/10/21  Bosanska   new inputs -> [1xnUE]struct UE_output (also output), [1x1]double uu, [1x1]double SNR
%                                      CQI evaluation for AWGN
%               2008/11/07  Bosanska   channel taps for PedA and PedB are shifted according to the
%                                      sampling freq LTE_params.Fs and rounded

ChanMod_output.genie.H_fft = zeros(LTE_params.Ntot, LTE_params.Nsub, ChanMod.nRX, ChanMod.nTX);
ChanMod_output.genie.H_fft_matrix = zeros(LTE_params.Ntot, LTE_params.Ntot, LTE_params.Nsub, ChanMod.nRX, ChanMod.nTX);

switch ChanMod.filtering
    case 'BlockFading'
        % generate ChanModel realization  
        switch ChanMod.type
            case {'AWGN'}
                %ChanMod_output.H = ones(ChanMod.nRX,ChanMod.nTX); % bad channel matrix, it has only rank one --> two codewords not possible
                H_temp = [  1,1,1,1;           
                            1,-1,-1,1;
                            1,-1,1,-1;
                            1,1,-1,-1];
%                 H_temp = exp(1i*2*pi*rand(4)); 
                ChanMod_output.H = H_temp(1:ChanMod.nRX,1:ChanMod.nTX); % now this channel matrix has rank = min(nRX,nTX)
                % Evaluation of the perfect channel knowledge in the frequency domain
                H_help = zeros(1,1,ChanMod.nRX, ChanMod.nTX);
                H_help(1,1,:,:) = ChanMod_output.H;
                ChanMod_output.genie.H_fft = repmat(H_help,[LTE_params.Ntot,LTE_params.Nsub,1,1]);

                % SINR and CQI evaluation for zero delay for the feedback channel
                SINR = LTE_common_calculate_SC_SNR( 'AWGN', SNR, LTE_params.Ntot );
            case {'flat Rayleigh'}
                if LTE_params.read_channel_from_trace
                    ChanMod_output.H = LTE_params.channel_matrix_trace.H_trace_normalized(:,:,subframe_i);
                else
                    ChanMod_output.H = 1/sqrt(2)*(randn(LTE_params.channel_param_RandStream,ChanMod.nRX,ChanMod.nTX)+1i*randn(LTE_params.channel_param_RandStream,ChanMod.nRX,ChanMod.nTX));
                end
                H_help(1,1,:,:) = ChanMod_output.H;
                ChanMod_output.genie.H_fft = repmat(H_help,[LTE_params.Ntot,LTE_params.Nsub,1,1]);
                %%%%%%%%%%%%%%%% this is not corrected, just to assign some
                %%%%%%%%%%%%%%%% values in order to get the simulator to
                %%%%%%%%%%%%%%%% run (Dasa)
                SINR = LTE_common_calculate_SC_SNR( 'flat Rayleigh', SNR, LTE_params.Ntot );
                %%%%%%%%%%%%%%%%
            case {'PedA', 'PedB', 'PedBcorr','VehA','VehB','TU','RA','HT','EPedA','EVehA','ETU'}
                
               % Generate faders
               switch ChanMod.time_correlation 
                   case 'correlated'
                       c = 299792458;
                       f = 2110e6;  % Frequency at which our system operates
                       v = UE.user_speed;  %speed at which we move
                       M = ChanMod.sin_num; %number of sin realizations
                       w_d = 2*pi*v*f/c;   % maximum radian Doppler frequency
                       time_i = (subframe_i-1)*LTE_params.Tsubframe;
                
                       % Yahong Rosa Zheng; Chengshan Xiao, "Simulation models with correct statistical properties for Rayleigh fading channels," IEEE Transactions on Communications, vol.51, no.6, pp. 920-928, June 2003
                       % URL:http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1209292&isnumber=27219
                       
                       % Zemen, T.; Mecklenbrauker, C.F., "Time-Variant Channel Estimation Using Discrete Prolate Spheroidal Sequences," IEEE Transactions on Signal Processing, vol.53, no.9, pp. 3597-3607, Sept. 2005
                       % URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1495893&isnumber=32146
                       
                       PI_mat = zeros(UE.nRX, ChanMod.nTX,ChanMod.interpolator.num_faders,M);
                       PI_mat(1,1,1,:) = (1:M)*2*pi;
                       PI_mat = repmat(PI_mat(1,1,1,:), [UE.nRX, ChanMod.nTX, ChanMod.interpolator.num_faders, 1]);
                       theta = zeros(UE.nRX, ChanMod.nTX,ChanMod.interpolator.num_faders,M);
                       theta(:,:,:,1) = UE.channel_coef_rosa_zheng.theta;
                       theta = repmat(theta(:,:,:,1), [1, 1, 1, M]);
                       alpha_n = (PI_mat - pi + theta) / (4*M);
                       
                       X_c = cos(UE.channel_coef_rosa_zheng.psi).*cos(w_d.*time_i.*cos(alpha_n) + UE.channel_coef_rosa_zheng.phi);
                       X_s = sin(UE.channel_coef_rosa_zheng.psi).*cos(w_d.*time_i.*cos(alpha_n) + UE.channel_coef_rosa_zheng.phi);
                       G = 2/sqrt(2*M) * sum(X_c + 1i*X_s,4);
                   case 'independent'
                       G = (randn(LTE_params.channel_param_RandStream,ChanMod.nRX,ChanMod.nTX,ChanMod.interpolator.num_faders) + 1i*randn(LTE_params.channel_param_RandStream,ChanMod.nRX,ChanMod.nTX,ChanMod.interpolator.num_faders)) /sqrt(2);
               end                 
                 
               % Generate channel matrix or read it from the trace
               if LTE_params.read_channel_from_trace
                   % Read channel matrix from trace
                   ChanMod_output.H = LTE_params.channel_matrix_trace.H_trace_normalized(:,:,:,subframe_i);
               else
                   % Generate and normalize H
                   ChanMod_output.H = ChanMod.interpolator.generateChannelMatrix(G,ChanMod.corrTX,ChanMod.corrRX);
                   ChanMod_output.H = ChanMod_output.H./ChanMod.normH;
               end
                
                % Evaluation of the perfect channel knowledge in the frequency domain
                for rr = 1:ChanMod.nRX
                    for tt = 1:ChanMod.nTX
                        spec = fft([squeeze(ChanMod_output.H(rr,tt,:)); zeros(LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,:))),1)]);                        
                        ChanMod_output.genie.H_fft(:,:,rr,tt) = repmat(spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]),1,LTE_params.Nsub); % remove DC carrier and zeros padded up to size of FFT
                    end
                end
                
                % CQI evaluation for zero delay for the feedback channel
                % !!!!!!!!!!!now works only for SISO, for MIMO it considers
                % the channel only for the first transmit and first receive antenna
                SINR = LTE_common_calculate_SC_SNR( 'frequency selective', SNR, LTE_params.Ntot,ChanMod_output.genie.H_fft );
            case 'winner_II'
                optargin = size(varargin,2);
                if optargin==2
                    user_channel_realization = varargin{1};
                    ChanMod_output.H = user_channel_realization;
                    % Evaluation of the perfect channel knowledge in the frequency domain
                    for rr = 1:ChanMod.nRX
                        for tt = 1:ChanMod.nTX
                            spec = fft([squeeze(ChanMod_output.H(rr,tt,:)); zeros(LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,:))),1)]);
                            ChanMod_output.genie.H_fft(:,:,rr,tt) = repmat(spec([LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 1:LTE_params.Ntot/2]),1,LTE_params.Nsub); % remove DC carrier and zeros padded up to size of FFT
                        end
                    end
                    
                    % CQI evaluation for zero delay for the feedback channel
                    % !!!!!!!!!!!now works only for SISO, for MIMO it considers
                    % the channel only for the first transmit and first receive antenna
                    SINR = LTE_common_calculate_SC_SNR( 'frequency selective', SNR, LTE_params.Ntot,ChanMod_output.genie.H_fft );
                elseif optargin==1
                    %do nothing
                else
                    error('something is wrong with the winner channel model');
                end
         end
    case 'FastFading'
        switch ChanMod.type
            case {'PedA', 'PedB', 'PedBcorr','VehA','VehB','TU','RA','HT','EPedA', 'EVehA','ETU','Rayleigh2'}
                % Yahong Rosa Zheng; Chengshan Xiao, "Simulation models with correct statistical properties for Rayleigh fading channels," Communications, IEEE Transactions on , vol.51, no.6, pp. 920-928, June 2003
                % URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1209292&isnumber=27219
                
                c = LTE_params.speed_of_light;
                f = LTE_params.carrier_freq;  % Frequency at which our system operates
                v = UE.user_speed;  %speed at which we move
                M = ChanMod.sin_num; %number of sin realizations
                w_d = 2*pi*v*f/c;   % maximum radian Doppler frequency
                ChanMod.interpolator.number_of_realizations = signal_length;
                number_of_realizations = signal_length;
                number_of_taps = ChanMod.interpolator.num_faders;
                X_c = zeros(UE.nRX, ChanMod.nTX,number_of_realizations,number_of_taps,M);
                X_s = zeros(UE.nRX, ChanMod.nTX,number_of_realizations,number_of_taps,M);
                psi_n = zeros(UE.nRX, ChanMod.nTX,number_of_realizations,number_of_taps,M);
                theta = zeros(UE.nRX, ChanMod.nTX,number_of_realizations,number_of_taps,M);
                phi = zeros(UE.nRX, ChanMod.nTX,number_of_realizations,number_of_taps,M);
                time_i = zeros(UE.nRX, ChanMod.nTX,number_of_realizations,number_of_taps,M);
                
                switch ChanMod.interpolation_method
                    case 'shift_to_nearest_neighbor'
                        ChanMod_output.H = zeros(UE.nRX, ChanMod.nTX, number_of_realizations, ChanMod.interpolator.num_faders);
                        time_i_help = unique(ChanMod.interpolator.tap_delays_samples)/ChanMod.interpolator.Fs;
                        time_i_help = kron(time_i_help, ones(number_of_realizations,1)) + (0:number_of_realizations-1).'*ones(1,length(time_i_help))*LTE_params.SamplingTime;
                        time_i(1,1,:,:,1) = time_i_help;
                        time_i = repmat(time_i(1,1,:,:,1),[UE.nRX, ChanMod.nTX, 1, 1, M]);
                        time_i = time_i + (subframe_i-1)*LTE_params.Tsubframe*ones(size(time_i));
                    case 'sinc_interpolation'
                        ChanMod_output.H = zeros(UE.nRX, ChanMod.nTX, number_of_realizations, length(ChanMod.interpolator.t));
                        time_i_help = ChanMod.interpolator.tap_delays;
                        time_i_help = kron(time_i_help, ones(number_of_realizations,1)) + (0:number_of_realizations-1).'*ones(1,length(time_i_help))*LTE_params.SamplingTime;
                        time_i(1,1,:,:,1) = time_i_help;
                        time_i = repmat(time_i(1,1,:,:,1),[UE.nRX, ChanMod.nTX, 1, 1, M]);
                        time_i = time_i + (subframe_i-1)*LTE_params.Tsubframe*ones(size(time_i));
                end
                                 
                switch ChanMod.time_correlation
                    case 'correlated'
                        psi_n(:,:,1, :, :)   = UE(uu).channel_coef_rosa_zheng.psi;
                        phi(:,:,1, :, :)     = UE(uu).channel_coef_rosa_zheng.phi;
                        theta(:,:,1, :, 1)   = UE(uu).channel_coef_rosa_zheng.theta;
                    case 'independent'
                        psi_n(:,:,1, :, :)   = rand(LTE_params.channel_param_RandStream,UE.nRX, ChanMod.nTX,number_of_taps, M);
                        phi(:,:,1, :, :)     = rand(LTE_params.channel_param_RandStream,UE.nRX, ChanMod.nTX,number_of_taps, M);
                        theta(:,:,1, :, 1)   = rand(LTE_params.channel_param_RandStream,UE.nRX, ChanMod.nTX,number_of_taps);
                        psi_n(:,:,1, :, :)   = (psi_n(:,:,1, :, :)*2 - 1) * pi;
                        phi(:,:,1, :, :)     = (phi(:,:,1, :, :)*2 - 1) * pi;
                        theta(:,:,1, :, 1)   = (theta(:,:,1, :, 1)*2 - 1) * pi;
                end
                
                cos_psi_n = cos(psi_n);
                sin_psi_n = sin(psi_n);
                cos_psi_n = repmat(cos_psi_n(:,:,1, :, :), [1,1,number_of_realizations, 1, 1]);
                sin_psi_n = repmat(sin_psi_n(:,:,1, :, :), [1,1,number_of_realizations, 1, 1]);
                psi_n = repmat(psi_n(:,:,1, :, :), [1,1,number_of_realizations, 1, 1]);
                phi = repmat(phi(:,:,1, :, :), [1,1,number_of_realizations, 1, 1]);
                theta = repmat(theta(:,:,1, :, 1), [1,1,number_of_realizations, 1, M]);
                
                PI_mat_minus_pi = zeros(UE.nRX, ChanMod.nTX,number_of_realizations,number_of_taps,M);
                PI_mat_minus_pi(1,1,1,1,:) = (1:M)*2*pi - pi;
                PI_mat_minus_pi = repmat(PI_mat_minus_pi(1,1,1,1,:), [UE.nRX, ChanMod.nTX, number_of_realizations, number_of_taps, 1]);
                alpha_n = (PI_mat_minus_pi + theta) / (4*M);
                
                cos_w_d_time_i = cos(w_d.*time_i.*cos(alpha_n) + phi);
                X_c = cos_psi_n.*cos_w_d_time_i;
                X_s = sin_psi_n.*cos_w_d_time_i;
                G = 2/sqrt(2*M) * sum(X_c + 1i*X_s,5);
                
                if strcmp(ChanMod.type,'Rayleigh2')
                    ChanMod_output.H = G;
                else
                    ChanMod_output.H = ChanMod.interpolator.generateChannelMatrix(G,ChanMod.corrTX,ChanMod.corrRX);
                end

                ChanMod_output.H = ChanMod_output.H./ChanMod.normH;
                 
                % NOTE: how can this be calculated for the Fast Fading case?
                SINR = 0;
                
                %% calculate FFT of channel matrix (necessary for feedback delay 0 simulations)
                for rr = 1:ChanMod.nRX
                    for tt = 1:ChanMod.nTX
                        taps_num = size(ChanMod_output.H,4);    %number of channel taps
                         if(strcmp(LTE_params.CyclicPrefix,'normal'))
                             start = 1;
                             for symbol_i = 1:LTE_params.Nsub
                                 if(symbol_i == 1 || symbol_i == 7)
                                     stop = start + LTE_params.NfftCP{1} + taps_num - 1 - 1;
                                     spec = fft([mean(squeeze(ChanMod_output.H(rr,tt,start:stop - taps_num + 1,:)),1), zeros(1, LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,1,:))))]);
                                     start = start + LTE_params.NfftCP{1};
                                 else
                                     stop = start + LTE_params.NfftCP{2} + taps_num - 1 - 1;
                                     spec = fft([mean(squeeze(ChanMod_output.H(rr,tt,start:stop - taps_num + 1,:)),1), zeros(1, LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,1,:))))]);
                                     start = start + LTE_params.NfftCP{2};
                                 end
                                 ChanMod_output.genie.H_fft(:,symbol_i,rr,tt) = spec(:,[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).';
                             end
                         else
                             start = 1;
                             for symbol_i = 1:LTE_params.Nsub
                                 stop = start + LTE_params.NfftCP + taps_num - 1 - 1;
                                 spec = fft([mean(squeeze(ChanMod_output.H(rr,tt,start:stop - taps_num + 1,:)),1), zeros(1, LTE_params.Nfft-length(squeeze(ChanMod_output.H(rr,tt,1,:))))]);
                                 start = start + LTE_params.NfftCP;
                                 ChanMod_output.genie.H_fft(:,symbol_i,rr,tt) = spec(:,[LTE_params.Nfft-LTE_params.Ntot/2+1:LTE_params.Nfft 2:LTE_params.Ntot/2+1]).';
                             end
                         end
                    end
                end
                                
            case 'winner_II'
                optargin = size(varargin,2);
                if optargin==2
                    user_channel_realization = varargin{1};
                    ChanMod_output.H = permute(user_channel_realization,[1 2 4 3]);
                    SINR = 0;
                elseif optargin==1
                    %do nothing
                else
                    error('something is wrong with the winner channel model');
                end
            otherwise
                error('Channel type not supported for fast fading');
        end
    otherwise
        error('chanMod.filtering type not supported');
end

% Store trace of the channel matrix for later reuse
if LTE_params.store_channel_trace
    LTE_params.channel_matrix_trace.store_H(ChanMod_output.H);
end
              
