function [ SINR ] = LTE_common_calculate_SC_SNR( chan_type, SNR, Ntot, varargin )
% Calculates the subcarrier SNR for the Block Fading cases.
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at. Original code From
% Dagmar Bosanska, dbosansk@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

sigma_n2 = 10^(-SNR/10);

switch chan_type
    case 'AWGN'
        SINR = 10*log10(1/sigma_n2) * ones(Ntot,1); % Flat, but still a per-subcarrier-SINR
    case 'flat Rayleigh'
        %%%%%%%%%%%%%%%% this is not corrected, just to assign some
        %%%%%%%%%%%%%%%% values in order to get the simulator to
        %%%%%%%%%%%%%%%% run (Dasa)
        SINR = 10*log10(1/sigma_n2) * ones(Ntot,1); % Flat, but still a per-subcarrier-SINR
    case 'frequency selective'
        genie.H_fft = varargin{1};
        % Set to work for SISO channels (2-dimensional matrix) and the
        % first symbol.
        SINR = 10*log10(abs(genie.H_fft(:,1)).^2/sigma_n2);
    otherwise
        error('Channel type not supported');
end
