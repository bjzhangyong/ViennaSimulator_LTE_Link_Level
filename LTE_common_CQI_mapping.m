function CQI = LTE_common_CQI_mapping(CQI_mapping_params,SINR)
% Map a SINR value to a floating-point CQI value (you will have to round it
% afterwards). Based on a linear interpolation taken from simulations. 
% [CQI] = LTE_common_CQI_mapping(SINR)
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2008 by INTHFT
% www.nt.tuwien.ac.at
%
% input:    SINR     ... [1 x 1]double - value or array of SINRs
% output:   CQI      ... [1 x 1]double or  [1 x length(SINR)]double - value or array of CQIs
%
% date of creation: 2008/10/21
% last changes

CQI = polyval(CQI_mapping_params.coeffs,SINR);
more_than_15 = (CQI>15);
less_than_0  = (CQI<0);
CQI(more_than_15) = 15;
CQI(less_than_0)  = 0;
