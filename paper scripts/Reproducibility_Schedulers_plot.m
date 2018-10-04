close all
%% LL
cd(LL_dir);

sum_TP2 = zeros(5,1);
j = 0;
confidence2 = zeros(5,2);

load('LL_BCQI.mat');
TP2 = zeros(length(simulation_results.UE_specific),5);
j = j+1;
for i = 1:length(simulation_results.UE_specific)
    TP2(i,j) = mean(sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
end
sum_TP2(j) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
confidence2(j,:) = bootci(length(sum(simulation_results.cell_specific.throughput_coded,3)),{@mean, sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6},'alpha',0.05); 

% load('Max_TP_LL.mat');
% j = j+1;
% for i = 1:length(simulation_results.UE_specific)
%     TP2(i,j) = mean(sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
% end
% sum_TP2(j) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
% confidence(j,:) = bootci(length(sum(simulation_results.cell_specific.throughput_coded,3)),{@mean, sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('LL_MM.mat');
j = j+1;
for i = 1:length(simulation_results.UE_specific)
    TP2(i,j) = mean(sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
end
sum_TP2(j) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
confidence2(j,:) = bootci(length(sum(simulation_results.cell_specific.throughput_coded,3)),{@mean, sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('LL_RR.mat');
j = j+1;
for i = 1:length(simulation_results.UE_specific)
    TP2(i,j) = mean(sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
end
sum_TP2(j) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
confidence2(j,:) = bootci(length(sum(simulation_results.cell_specific.throughput_coded,3)),{@mean, sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('LL_PF.mat');
j = j+1;
for i = 1:length(simulation_results.UE_specific)
    TP2(i,j) = mean(sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
end
sum_TP2(j) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
confidence2(j,:) = bootci(length(sum(simulation_results.cell_specific.throughput_coded,3)),{@mean, sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('LL_RF.mat');
j = j+1;
for i = 1:length(simulation_results.UE_specific)
    TP2(i,j) = mean(sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6;
end
sum_TP2(j) = mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6;
confidence2(j,:) = bootci(length(sum(simulation_results.cell_specific.throughput_coded,3)),{@mean, sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6},'alpha',0.05); 



[SNR_vec_sort,indices] = sort(SNR_vec);

figure(4);
h3 = bar3(TP2(indices,:));
% set(h3,'facealpha',1)
set(gca,'Fontsize',24)
ylim([0,21])
ylabel('User index')
% ylab = get(gca,'Ylabel');
% set(ylab,'Position',[-62.4548  -97.5486    6.5281]);
zlabel('Throughput [Mbit/s]');
l= legend('Best CQI','Max. min.','Round robin','Prop. fair','Res. fair');
% set(l,'Position',[0.2583    0.5703    0.1437    0.2106]);
% set(gca,'view',[-149.5000    8.0000]);
for jj = 1:j
    J2(jj) = sum(TP2(:,jj))^2/(length(simulation_results.UE_specific)*sum(TP2(:,jj).^2));
end
title('LL throughputs');

%% SL
cd(SL_dir);

sum_TP = zeros(5,1);
confidence = zeros(5,2);
j = 0;

load('SL_BCQI.mat');
TP = zeros(20,5);
% SNR = zeros(20,5);
% SINR = zeros(20,5);
j = j+1;
for i = 1:20
    TP(i,j) = sum(simulation_traces.UE_traces(i).TB_size(1,:).*uint32(simulation_traces.UE_traces(i).ACK(1,:)))/size(simulation_traces.UE_traces(i).TB_size,2)*10^-3;
%     SNR(i,j) = simulation_traces.UE_traces(i).SNR(1);
%     SINR(i,j) = 10*log10(mean(mean(squeeze(10.^(simulation_traces.UE_traces(i).SINR(1,:,:)/10)))));
end
sum_TP(j) = sum(TP(:,j));
sector_TP = simulation_traces.eNodeB_tx_traces.sector_traces(1).acknowledged_data(1,:);
confidence(j,:) = bootci(length(sector_TP),{@mean, sector_TP/LTE_params.Tsubframe/1e6},'alpha',0.05); 

% load('Max_TP.mat');
% j = j+1;
% for i = 1:20
%     TP(i,j) = sum(simulation_traces.UE_traces(i).TB_size(1,:))/size(simulation_traces.UE_traces(i).TB_size,2)*10^-3;
%     SNR(i,j) = simulation_traces.UE_traces(i).SNR(1);
%     SINR(i,j) = 10*log10(mean(mean(squeeze(10.^(simulation_traces.UE_traces(i).SINR(1,:,:)/10)))));
% end
% sum_TP(j) = sum(TP(:,j));
% sector_TP = simulation_traces.eNodeB_tx_traces.sector_traces(j).acknowledged_data(1,:);
% confidence(j,:) = bootci(length(sector_TP),{@mean, sector_TP/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('SL_MM.mat');
j = j+1;
for i = 1:20
     TP(i,j) = sum(simulation_traces.UE_traces(i).TB_size(1,:).*uint32(simulation_traces.UE_traces(i).ACK(1,:)))/size(simulation_traces.UE_traces(i).TB_size,2)*10^-3;
end
sum_TP(j) = sum(TP(:,j));
sector_TP = simulation_traces.eNodeB_tx_traces.sector_traces(1).acknowledged_data(1,:);
confidence(j,:) = bootci(length(sector_TP),{@mean, sector_TP/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('SL_RR.mat');
j = j+1;
for i = 1:20
    TP(i,j) = sum(simulation_traces.UE_traces(i).TB_size(1,:).*uint32(simulation_traces.UE_traces(i).ACK(1,:)))/size(simulation_traces.UE_traces(i).TB_size,2)*10^-3;
%     SNR(i,j) = simulation_traces.UE_traces(i).SNR(1);
%     SINR(i,j) = 10*log10(mean(mean(squeeze(10.^(simulation_traces.UE_traces(i).SINR(1,:,:)/10)))));
end
sum_TP(j) = sum(TP(:,j));
sector_TP = simulation_traces.eNodeB_tx_traces.sector_traces(1).acknowledged_data(1,:);
confidence(j,:) = bootci(length(sector_TP),{@mean, sector_TP/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('SL_PF.mat');
j = j+1;
for i = 1:20
    TP(i,j) = sum(simulation_traces.UE_traces(i).TB_size(1,:).*uint32(simulation_traces.UE_traces(i).ACK(1,:)))/size(simulation_traces.UE_traces(i).TB_size,2)*10^-3;
%     SNR(i,j) = simulation_traces.UE_traces(i).SNR(1);
%     SINR(i,j) = 10*log10(mean(mean(squeeze(10.^(simulation_traces.UE_traces(i).SINR(1,:,:)/10)))));
end
sum_TP(j) = sum(TP(:,j));
sector_TP = simulation_traces.eNodeB_tx_traces.sector_traces(1).acknowledged_data(1,:);
confidence(j,:) = bootci(length(sector_TP),{@mean, sector_TP/LTE_params.Tsubframe/1e6},'alpha',0.05); 

load('SL_RF.mat');
j = j+1;
for i = 1:20
    TP(i,j) = sum(simulation_traces.UE_traces(i).TB_size(1,:).*uint32(simulation_traces.UE_traces(i).ACK(1,:)))/size(simulation_traces.UE_traces(i).TB_size,2)*10^-3;
%     SNR(i,j) = simulation_traces.UE_traces(i).SNR(1);
%     SINR(i,j) = 10*log10(mean(mean(squeeze(10.^(simulation_traces.UE_traces(i).SINR(1,:,:)/10)))));
end
sum_TP(j) = sum(TP(:,j));
sector_TP = simulation_traces.eNodeB_tx_traces.sector_traces(1).acknowledged_data(1,:);
confidence(j,:) = bootci(length(sector_TP),{@mean, sector_TP/LTE_params.Tsubframe/1e6},'alpha',0.05); 

% [SNR_vec_sort,indices] = sort(SINR(:,1));

figure(3);
h3 = bar3(TP(indices,:));
set(gca,'Fontsize',24)
ylim([0,21])
ylabel('User index')
% ylab = get(gca,'Ylabel');
% set(ylab,'Position',[-62.4548  -97.5486    6.5281]);
zlabel('Throughput [Mbit/s]');
l= legend('Best CQI','Max. min.','Round robin','Prop. fair','Res. fair');
% set(l,'Position',[0.2583    0.5703    0.1437    0.2106]);
% set(gca,'view',[-149.5000    8.0000]);
title('SL throughputs');

for jj = 1:j
    J(jj) = sum(TP(:,jj))^2/(20*sum(TP(:,jj).^2));
end

figure(1)
h = bar([J',J2'],1);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
set(gca,'Fontsize',24)
legend('System level','Link level','Location','Northwest')
ylabel('Jain´s fairness index');
set(gca,'XTicklabel',{'Best CQI','Max. min.','Round robin','Prop. fair','Res. fair'});
grid on 

figure(2)
h = bar([sum_TP,sum_TP2],1);
set(h(1),'facecolor','red')
set(h(2),'facecolor','blue')
set(gca,'Fontsize',24)
legend('System level','Link level','Location','Northwest')
ylabel('Throughput [Mbit/s]');
set(gca,'XTicklabel',{'Best CQI','Max. min.','Round robin','Prop. fair','Res. fair'});
grid on
hold on

plot_color = 'green';
for index = 1:length(confidence)
    plot([index-0.125 index-0.125], [confidence(index,1) confidence(index,2)],'Color',plot_color,'Markersize',10,'Linewidth',2.5);
    plot(index-0.125+[-0.1 0.1], [confidence(index,1) confidence(index,1)],'Color',plot_color,'Markersize',10,'Linewidth',2.5);
    plot(index-0.125+[-0.1 0.1], [confidence(index,2) confidence(index,2)],'Color',plot_color,'Markersize',10,'Linewidth',2.5);
    
    plot([index+0.125 index+0.125], [confidence2(index,1) confidence2(index,2)],'Color',plot_color,'Markersize',10,'Linewidth',2.5);
    plot(index+0.125+[-0.1 0.1], [confidence2(index,1) confidence2(index,1)],'Color',plot_color,'Markersize',10,'Linewidth',2.5);
    plot(index+0.125+[-0.1 0.1], [confidence2(index,2) confidence2(index,2)],'Color',plot_color,'Markersize',10,'Linewidth',2.5);
end

% figure(5)
% bar(SNR_vec_sort)
% set(gca,'Fontsize',24)
% ylabel('SNR [dB]');
% xlabel('User index');
% grid on
% xlim([0,21])

% save('Scheduler_data','TP','TP2','J','J2','confidence','confidence2','SNR_vec_sort','indices');