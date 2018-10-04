function EESMs = LTE_common_EESM( SINRs, beta,mode )
% Calculate EESM from the input column vector/matrix. EESM is calculated
% for each column of the N-dimensional input matrix. Input can either be
% linear or logarithmic (dB). In the dB case, the value is converted to
% linear before the EESM averaging and then converted again to dB before
% returning the results.
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at
%
% For more information, you can see:
% Blankenship, Y.W.; Sartori, P.J.; Classon, B.K.; Desai, V.; Baum, K.L., "Link error prediction methods for multicarrier systems," Vehicular Technology Conference, 2004. VTC2004-Fall. 2004 IEEE 60th , vol.6, no., pp. 4175-4179 Vol. 6, 26-29 Sept. 2004
% URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1404865&isnumber=30415

switch mode
    case 'dB'
        in_dB = true;
        SINRs = 10.^(SINRs/10);
    case 'lin'
        in_dB = false;
    otherwise
        error('Only "lin" or "dB" accepted"');
end

length_SINRs = prod(size(SINRs));
step = size(SINRs,1);
EESMs = zeros(length_SINRs/step,1);

j_ = 1;
for i_=1:step:length_SINRs
    SINRS_to_average = SINRs(i_:i_+step-1);
    EESMs(j_) = -beta*log(mean(exp(-SINRS_to_average/beta)));
    j_ = j_+1;
end

if in_dB
    EESMs = 10*log10(EESMs);
end