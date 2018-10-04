clear TP
figure(2)
clf
clc
clear y
N_UE = length(simulation_results.UE_specific);
av_wind = 200;
% a = [1,-(1-1/av_wind)];
% b = 1/av_wind;
a = 1;
b = 1/av_wind *ones(av_wind,1); % moving average filter
color = [1,0,0;0,0,1;0,1,0;1,0,1;0,1,1];
for i = 1:N_UE
    y(i,:) = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded(1:end,1,:),3));
    TP(i) = y(i,end)/LTE_params.Tsubframe/1e6;
%     TP(i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
    h(i) = plot(y(i,:)/LTE_params.Tsubframe/1e6,'Color',color(i,:),'Linewidth',2.5);
%     color = circshift(color,[0,-1]);
    hold on
    grid on
end
set(gca,'Fontsize',30)
% set(gca,'Ylim',[0,0.6],'Ytick',[0:0.1:0.6]);
set(gca,'Ylim',[0,0.7],'Ytick',[0:0.1:0.7]);
ylabel('Throughput [Mbit/s]');
pos = get(gca,'Position');
ax2 = axes('position',pos);
J = sum(y,1).^2./(N_UE*sum(y.^2,1));
h(i+1) = plot(ax2,J,'k','Linewidth',2.5);
set(ax2,'color','none','YAxisLocation','right','Ticklength',[0.02,0.025]);
set(gca,'Fontsize',30)
% set(gca,'Ylim',[0.5,0.85],'Ytick',[0.5:0.05:0.85]);
xlabel('Subframe number')
ylabel('Jain\primes fairness index');
legend(h,'Throughput User 1','Throughput User 2','Throughput User 3','Throughput User 4','Throughput User 5','Fairness','Location','Southeast')

for i = 1:N_UE
   LT_TP(i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
end
LT_J = sum(LT_TP)^2/(N_UE*sum(LT_TP.^2))
% J = sum(TP)^2/(N_UE*sum(TP.^2));
% J
% TP
% sum(TP) 