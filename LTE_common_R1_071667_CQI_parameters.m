function CQI_parameters = LTE_common_R1_071667_CQI_parameters
% It returns the MCSs used in the R1-07196 RAN document. Used for
% result validation purposes
% Author: Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

rates = [ 1/9 1/6 0.21 1/4 1/3 0.42 1/2 0.58 2/3 0.73 0.43 0.46 1/2 0.54 0.58 0.61 2/3 0.73 4/5 0.58 0.62 2/3 0.70 0.74 4/5 0.85 0.90];
modulation_orders = [2*ones(1,10) 4*ones(1,9) 6*ones(1,8)];

for i_=1:length(rates)
    CQI_parameters(i_).CQI = 100+i_;
    switch modulation_orders(i_)
        case 2
            CQI_parameters(i_).modulation = 'QPSK';
        case 4
            CQI_parameters(i_).modulation = '16QAM';
        case 6
            CQI_parameters(i_).modulation = '64QAM';
    end
    CQI_parameters(i_).modulation_order = modulation_orders(i_);
    CQI_parameters(i_).coding_rate_x_1024 = rates(i_)*1024;
    CQI_parameters(i_).efficiency = modulation_orders(i_)*rates(i_);
end
