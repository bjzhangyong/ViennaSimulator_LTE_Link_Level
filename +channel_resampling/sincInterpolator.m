classdef sincInterpolator < handle
% Implementation a sinc interpolator. The acausal part of the generated impulse
% response can also be used (configure include_acausal_part=true for that).
% (c) Josep Colom Ikuno, INTHFT, 2009

properties
    tap_delays              % Tap delays
    relative_power_dB       % Relative powers (dB)
    
    PDPFreq                 % Sampling frequency of the Power Density Function
    norm_h                  % h normalization factor
    t                       % range for which h is generated (from -range to +range)
    
    number_of_taps          % Number of taps of this model
    num_faders              % Number of faders needed. In this case, this is equal to the number of taps (number of sincs)
    
    precomputed_sincs       % pregenerated sincs used for interpolation
    
    nTX                     % Number of TX antennas
    nRX                     % Number of RX antennas
    
    number_of_realizations  % Number of channel realizations, mainly for fast fading case
end

methods
    % Class constructor
    %   - fg: bandwidth of your channel model
    %   - the_tap_delays_in_s: where your taps are
    %   - the_relative_powers_in_dB: power of the taps (in dB)
    %   - include_acausal_part: whether you want to include the acausal
    %       part in your generated sincs.
    %   - num_samples: how many samples long will the sinc-interpolated
    %       channel impulse response be (the positive part). If you set h
    %       to be acausal, then h will be 2*num_samples-1 samples long.
    %   - fader_source: how the random faders will be generated
    function obj = sincInterpolator(...
            f_g,...
            the_tap_delays_in_s,...
            the_relative_powers_in_dB,...
            include_acausal_part,...
            num_samples,...
            nTX,nRX)
        
        obj.tap_delays = the_tap_delays_in_s;
        obj.relative_power_dB   = the_relative_powers_in_dB;
        obj.PDPFreq = f_g;
        obj.nTX = nTX;
        obj.nRX = nRX;

        if length(obj.tap_delays)~=length(obj.relative_power_dB)
            error('Every tap delay power gain must be specified.');
        end

        number_of_taps = length(the_tap_delays_in_s);
        obj.number_of_taps = number_of_taps;
        obj.num_faders = number_of_taps;

        % number of samples (positive side)
        if include_acausal_part
            t = (-num_samples:num_samples) / f_g;
        else
            t = (0:num_samples) / f_g;
        end

        % Allocate sincs
        obj.precomputed_sincs =  zeros(number_of_taps,length(t));

        % The sinc function is defined in Matlab as
        %   sinc(t) = sin(pi*t)/(pi*t)
        for i_=1:number_of_taps
            obj.precomputed_sincs(i_,:) = sinc(2*f_g*(t-the_tap_delays_in_s(i_))) .* sqrt(10.^(0.1*obj.relative_power_dB(i_)));
        end

        obj.norm_h = sum(sum(obj.precomputed_sincs).^2);
        obj.t = t;
    end
    
    % Multiply the specified faders with the sincs and add them to get a
    % channel impulse response
    %   -If acausal=true it will give you h including the acausal part
    %   -If acausal=false h will NOT include the acausal part
    function [h] = generateChannelMatrix(obj,G,corrTX,corrRX)
        
        if ndims(G)==3  %block fading case
            
            h_temp = zeros(obj.nRX,obj.nTX,obj.number_of_taps,length(obj.t));
            h = zeros(obj.nRX,obj.nTX,length(obj.t));

            % Need to add the correlation matrices, which right now are not the
            % correct size
            for tap_idx=1:obj.number_of_taps
                for rx_ant = 1:size(G,1)
                    for tx_ant = 1:size(G,2)
                        h_temp(rx_ant,tx_ant,tap_idx,:) = G(rx_ant,tx_ant,tap_idx) * obj.precomputed_sincs(tap_idx,:);
                    end
                end
            end
            h(:,:,:) = squeeze(sum(h_temp,3));  % Sum all over the sincs to obtain the final matrix

            h = h / obj.norm_h; % NOTE: please check normalization here
            
        elseif ndims(G)==4  %fast fading case
            
            h_temp = zeros(obj.nRX,obj.nTX,obj.number_of_realizations,obj.number_of_taps,length(obj.t));
            h = zeros(obj.nRX,obj.nTX,obj.number_of_realizations,length(obj.t));

            % Need to add the correlation matrices, which right now are not the
            % correct size
            for tap_idx=1:obj.number_of_taps
                for rx_ant = 1:size(G,1)
                    for tx_ant = 1:size(G,2)
                        h_temp(rx_ant,tx_ant,:,tap_idx,:) = squeeze(G(rx_ant,tx_ant,:,tap_idx)) * squeeze(obj.precomputed_sincs(tap_idx,:));
                    end
                end
            end
            h(:,:,:,:) = squeeze(sum(h_temp,4));  % Sum all over the sincs to obtain the final matrix

            h = h / obj.norm_h; % NOTE: please check normalization here
        end
    end
    
    % Generate h matrix without applying gains (only the sincs)
    function [h] = generateChannelMatrix_no_tap_gain(obj)
        h = 0;

        for i_=1:obj.number_of_taps
            h = h + obj.precomputed_sincs(i_,:);
        end
        
        h = h / obj.norm_h;
    end
    
    % Plot pre-generated sincs
    function plot_sincs(obj)
        figure1 = figure;
        box('on');
        grid on;
        hold('all');
        for i_=1:size(obj.precomputed_sincs,1)
            stem(obj.t*1e9,obj.precomputed_sincs(i_,:),'.');
        end
        xlabel({'time (ns)'});
        ylabel({'x(t)'});
        title(['pregenerated sincs, ' num2str(obj.PDPFreq/1e6,'%3.2f') ' MHz'])
    end
    
    % Plot a given generated h
    function plot_h(obj,h)
        figure1 = figure;
        axes('Parent',figure1);
        box('on');
        grid on;
        hold('all');
        stem(obj.t*1e9,abs(h),'.');
        ylim([0 max(abs(h))]);
        xlabel({'time (ns)'});
        ylabel({'|h(t)| (dB)'});
        title(['channel impulse response (h), ' num2str(obj.PDPFreq/1e6,'%3.2f') ' MHz']);
    end
    
    % Plot a given generated H
    function plot_H(obj,h)
        figure1 = figure;
        H = fftshift(fft(h));
        frequency_scale = linspace(-obj.PDPFreq/1e6,obj.PDPFreq/1e6,length(H));
        box('on');
        grid on;
        hold('all');
        plot(frequency_scale,10*log10(abs(H)));
        xlabel({'f (MHz)'});
        ylabel({'|H(f)|'});
        title(['channel impulse response (H), ' num2str(obj.PDPFreq/1e6,'%3.2f') ' MHz (from -\pi to \pi)'])
    end
    
end
end
