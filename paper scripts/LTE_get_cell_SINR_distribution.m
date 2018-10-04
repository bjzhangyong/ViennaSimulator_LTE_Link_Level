% Execute this script during a breakpoint positioned in the last line of
% the LTE SL simulator funcition LTE_common_calculate_cell_capacity

eNodeB_id = 11;
eNodeB_sector_area  = networkPathlossMap.sector_assignment_no_shadowing(:,:,2)==eNodeB_id; % Take all sector (better statistics)
eNodeB_sector_SINRs_cut = max_SINR_dB_all(eNodeB_sector_area);

% Reduce the number of points to a more reasonable number. ~200
N_SINRs = length(eNodeB_sector_SINRs_cut);
delmation = floor(N_SINRs/200);
eNodeB_sector_SINRs_cut_less = eNodeB_sector_SINRs_cut(1:delmation:end);

figure;
imagesc(eNodeB_sector_area);
set(gca,'YDir','normal');
title(sprintf('Area pertaining eNodeB %d, sector %d',eNodeB_id,sector_id));

figure;
hist(eNodeB_sector_SINRs_cut(:),100);
title('SINR hist');

figure;
hold all;
ecdf(eNodeB_sector_SINRs_cut);
ecdf(eNodeB_sector_SINRs_cut_less);
grid on;
legend('show','Location','Best',{'no delmation' sprintf('1/%d delmation',delmation)});

save('HARQ_validation_SINRs.mat','eNodeB_sector_SINRs_cut','eNodeB_sector_SINRs_cut_less');