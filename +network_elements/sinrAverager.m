classdef sinrAverager < handle
    % Defines the abstract classes needed by a class that implements a SINR
    % averaging method
    % (c) Josep Colom Ikuno, INTHFT, 2008
    % adopted for Link Level by Stefan Schwarz 2010

   properties
   end

   methods (Abstract)
       % INPUT VALUES IN LINEAR
       % OUTPUT VALUE IN dB
       effective_SINR = average(SINR_vector,varargin)
   end
end 
