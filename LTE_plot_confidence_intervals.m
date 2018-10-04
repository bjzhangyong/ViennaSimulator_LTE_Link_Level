function LTE_plot_confidence_intervals(SNR_vector, plot_color, figure_name,function_type,varargin)
% Plots all the meaningful data from the results file
% Michal Simko, msimko@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at


optargin = size(varargin,2);
if optargin==1
    data1 = varargin{1};
elseif optargin==2
    data1 = varargin{1};
    data2 = varargin{2};
else
    error('Wrong number of input variables');
end

figure(figure_name)
if strcmp('mean',function_type)
    for i=1:size(data1,2)
        temp = bootci(2000,{@mean, data1(:,i)},'alpha',0.05); 
        T_min(i) = temp(1);
        T_max(i) = temp(2);
    end
elseif strcmp('ber',function_type)
    ber_func = @(bit_error,nr_of_bits) sum(bit_error)/sum(nr_of_bits);
    for i=1:size(data1,2)
        temp = bootci(2000,{ber_func, data1(:,i),data2(:,i)},'alpha',0.05);
        T_min(i) = temp(1);
        if T_min(i) == 0
            T_min(i) = eps;
        end
        T_max(i) = temp(2);
    end
else
    error('wrong function type');
end
for index = 1:length(SNR_vector)
    plot([SNR_vector(index) SNR_vector(index)], [T_min(index) T_max(index)],'Color',plot_color,'Markersize',5);
    plot(SNR_vector(index)+[-0.2 0.2], [T_min(index) T_min(index)],'Color',plot_color,'Markersize',5);
    plot(SNR_vector(index)+[-0.2 0.2], [T_max(index) T_max(index)],'Color',plot_color,'Markersize',5);
end