classdef channelMatrixTrace < handle
% Class that store a trace of channel matrices, as well as the seed
% used for the noise generation.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at
    
    properties
        noise_seed          % Seed used for the noise
        H_trace_normalized  % Normalized channel matrices
        last_insertion = 0;
        
        % Things that will be checked when the trace is load to assert its validity
        bandwidth           % System Bandwidth
        type                % Type of channel
        nTx                 % Number of Transmit Antennas
        nRx                 % Number of Receive Antennas
        TTI_length          % Length of the trace in TTIs
        Nfft                % Number of points of the FFT
    end
    
    methods
        % Store a channel matrix
        function store_H(obj,H_normalized)
            if isempty(obj.H_trace_normalized)
                H_trace_normalized = zeros([ size(H_normalized) obj.TTI_length ]);
            end
            switch ndims(H_normalized)
                case 2 % AWGN and flat Rayleigh cases
                    obj.H_trace_normalized(:,:,obj.last_insertion+1)   = H_normalized;
                case 3 % Block fading case
                    obj.H_trace_normalized(:,:,:,obj.last_insertion+1)   = H_normalized;
                case 4 % Fast fading case
                    obj.H_trace_normalized(:,:,:,:,obj.last_insertion+1) = H_normalized;
                otherwise
                    error('Channel matrix dimensions do not match');
            end
            obj.last_insertion = obj.last_insertion + 1;
        end
        
        % Plot channel trace
        function plot_trace(obj)
            H_trace_normalized_fft = fftshift(fft(obj.H_trace_normalized,obj.Nfft,3),3);
            for nTx = 1:obj.nTx
                for nRx = 1:obj.nRx
                    data_to_plot = squeeze(abs(H_trace_normalized_fft(nRx,nTx,:,:))).^2;
                    
                    the_figure = figure;
                    the_axes = axes('Parent',the_figure);
                    the_title = sprintf('|H|^2 trace. TX antenna %d, RX antenna %d, %3.2f MHz bandwidth',nTx,nRx,obj.bandwidth/1e6);
                    
                    surf(the_axes,data_to_plot,'LineStyle','none');
                    title(the_axes,the_title);
                    xlim(the_axes,[1 obj.TTI_length]);
                    xlabel(the_axes,'TTI');
                    ylim(the_axes,[1 obj.Nfft]);
                    ylabel(the_axes,'FFT point');
                    zlabel('|H|^2');
                    %view(the_axes,[-21.5 58]);
                    grid(the_axes,'on');
                    %set(the_axes,'ZScale','log');
                    set(the_axes,'ZScale','linear');
                end
            end
        end
    end
    
end

