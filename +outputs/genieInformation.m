classdef genieInformation < handle
% Genie information. Including transmitted bits and so on.
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

   properties
       data_bits
       bits_to_turboencode
       sent_bits
       y_tx_assembled      % Contains the modulated sent symbols before the padding and IFFT
       v                   % added noise at RX-antennas (preFFT, y = H*x + v) with respect to Nfft and Ntot
       n                   % noise after FFT at the receiver (postFFT), where v = sqrt(LTE_params.Nfft/LTE_params.Ntot) * n
       layer_x             % contains 
   end

   methods
       function clear(obj)
           obj.data_bits = [];
           obj.bits_to_turboencode = [];
           obj.sent_bits = [];
       end
   end
end 
