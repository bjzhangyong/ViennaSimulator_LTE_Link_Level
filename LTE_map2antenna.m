function y_tx = LTE_map2antenna(y_tx_temp,mapping)
% This function maps antenna ports to antennas
% now there is only a dummy mapping
% Author: Stefan Schwarz, sschwarz@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at
y_tx(:,mapping) = y_tx_temp;