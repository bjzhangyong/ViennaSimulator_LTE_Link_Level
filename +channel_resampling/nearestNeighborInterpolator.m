classdef nearestNeighborInterpolator < handle
% Implementation a nearest neighbor interpolator.
% (c) Josep Colom Ikuno, INTHFT, 2009

   properties
       tap_delays               % Tap delays
       tap_delays_samples       % Tap delays in samples for this sampling frequency
       relative_power_dB        % Relative powers (dB)
       NTap                     % Number of taps
       Fs                       % New sampling frequency
       nTX                      % Number of TX antennas
       nRX                      % Number of RX antennas
       num_faders               % number of faders that we will need. As taps may merge, this may not be the number of taps the PDP has
       number_of_realizations   % Number of channel realizations, mainly for fast fading case
   end

   methods
       % Class constructor
       function obj = nearestNeighborInterpolator(...
               Fs,...
               NTap,...
               the_tap_delays_in_s,...
               the_relative_powers_in_dB,...
               nTX,nRX)
               
           obj.nTX = nTX;
           obj.nRX = nRX;
           obj.NTap = NTap;
           obj.Fs = Fs;
           obj.tap_delays = the_tap_delays_in_s;
           obj.tap_delays_samples = round(the_tap_delays_in_s*Fs);
           obj.relative_power_dB = the_relative_powers_in_dB;
           obj.num_faders = length(unique(obj.tap_delays_samples)); % some taps merge by a specific sampling freq
       end
       
       function [h] = generateChannelMatrix(obj,G,corrTX,corrRX)
           
           if ndims(G)==3  %block fading case
               
               h = zeros(obj.nRX,obj.nTX,obj.tap_delays_samples(end)+1);
               distinct_tap_delays_samples = unique(obj.tap_delays_samples);

               for tap_i = 1:obj.num_faders
                   equal_tap_delays_samples = distinct_tap_delays_samples(tap_i)==obj.tap_delays_samples;
                   h(:,:,distinct_tap_delays_samples(tap_i)+1) = ...
                       sqrt(sum(10.^(obj.relative_power_dB(equal_tap_delays_samples)./10))).*...
                       (sqrtm(squeeze(corrRX(tap_i,:,:)))*G(:,:,tap_i)*...
                       (sqrtm(squeeze(corrTX(tap_i,:,:)))).');
               end
               
           elseif ndims(G)==4  %fast fading case
               % NOTE: antenna correlation not implemented
               
               h = zeros(obj.nRX,obj.nTX,obj.number_of_realizations,obj.tap_delays_samples(end)+1);
               distinct_tap_delays_samples = unique(obj.tap_delays_samples);

               for tap_i = 1:obj.num_faders
                   equal_tap_delays_samples = distinct_tap_delays_samples(tap_i)==obj.tap_delays_samples;

                   h(:,:,:,distinct_tap_delays_samples(tap_i)+1) = ...
                       sqrt(sum(10.^(obj.relative_power_dB(equal_tap_delays_samples)./10))).*...
                       G(:,:,:,tap_i);
               end
               
           end
       end
   end 
end
