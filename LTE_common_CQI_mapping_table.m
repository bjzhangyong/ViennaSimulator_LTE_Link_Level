function CQI = LTE_common_CQI_mapping_table(CQI_mapping_params,SINRs,CQIs)
% The function takes as input argument a number of SINRs corresponding
% to the input set of CQIs and delivers as output the highest possible CQI 
% with BLER <= 0.1

% CQI = ones(size(SINR));
% for i1 = 1:size(SINR,1)
%     for i2 = 1:size(SINR,2)
%         temp = find(SINR(i1,i2) >= CQI_mapping_params.table,1,'last');
%         if isempty(temp)
%             temp = 1;
%         end
%         CQI(i1,i2) = temp;      
%     end
% end

temp = zeros(size(CQIs));
temp(CQI_mapping_params.table(CQIs) <= SINRs) = 1;
CQI = find(temp,1,'last')-1;
CQI(CQI == 0) = 20;

% if isempty(CQI)
%     CQI = 0;
% end
