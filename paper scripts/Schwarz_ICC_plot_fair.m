close all
cd ..
fairn = 0.5:0.025:0.975;
J = zeros(length(fairn),1);
TP_UE = zeros(length(fairn),2);
for fairi2 = 1:length(fairn)
    load(['SOCP_' num2str(fairi2) '.mat']);
    LTE_params.scheduler.av_window = 500;
    a = [1,-(1-1/LTE_params.scheduler.av_window)];
    b = 1/LTE_params.scheduler.av_window;
    N_UE = length(simulation_results.UE_specific);
    for i = 1:N_UE
%         y = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded,3));
%         TP_UE(fairi2,i) = y(end)/LTE_params.Tsubframe/1e6;
        TP_UE(fairi2,i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
    end
    J(fairi2) = sum(TP_UE(fairi2,:))^2/(N_UE*sum(TP_UE(fairi2,:).^2));
%     plot(J(fairi2),sum(TP(fairi2,:)),'gx','Markersize',14,'linewidth',2.5);
end
plot(J,sum(TP_UE,2),'bx-','Markersize',18,'linewidth',2.5);
hold on
grid on
set(gca,'Fontsize',30);
xlabel('Jain`s fairness index');
ylabel('Throughput [Mbit/s]');

fairn = [0,0.2,0.5,0.9,1.5,2,4,10,20,1000];
J = zeros(length(fairn),1);
TP_UE = zeros(length(fairn),2);
for fairi2 = 1:length(fairn)
    load(['Alpha_' num2str(fairi2) '.mat']);
    LTE_params.scheduler.av_window = 500;
    a = [1,-(1-1/LTE_params.scheduler.av_window)];
    b = 1/LTE_params.scheduler.av_window;
    N_UE = length(simulation_results.UE_specific);
    for i = 1:N_UE
%         clear y;
%         y = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded,3));
%         TP_UE(fairi2,i) = y(end)/LTE_params.Tsubframe/1e6;
        TP_UE(fairi2,i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
    end
    J(fairi2) = sum(TP_UE(fairi2,:))^2/(N_UE*sum(TP_UE(fairi2,:).^2));
%     plot(J(fairi2),sum(TP(fairi2,:)),'r+','Markersize',14,'linewidth',2.5);
%     hold on
%     grid on
end
sum_TP = sum(TP_UE,2);
plot(J,sum_TP,'r+-','Markersize',14,'linewidth',2.5);

load('MM.mat');
LTE_params.scheduler.av_window = 500;
    a = [1,-(1-1/LTE_params.scheduler.av_window)];
    b = 1/LTE_params.scheduler.av_window;
for i = 1:N_UE
%     y = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded,3));
%     TP_mm(i) = y(end)/LTE_params.Tsubframe/1e6;
    TP_mm(i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
end
J_mm = sum(TP_mm)^2/(N_UE*sum(TP_mm.^2));
plot(J_mm,sum(TP_mm),'ko','Markersize',16,'linewidth',2.5);

load('BCQI.mat');
LTE_params.scheduler.av_window = 500;
    a = [1,-(1-1/LTE_params.scheduler.av_window)];
    b = 1/LTE_params.scheduler.av_window;
for i = 1:N_UE
%     y = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded(1:5000,:,:),3));
%     TP_bc(i) = y(end)/LTE_params.Tsubframe/1e6;
    TP_bc(i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
end
J_bc = sum(TP_bc)^2/(N_UE*sum(TP_bc.^2));
plot(J_bc,sum(TP_bc),'ks','Markersize',16,'linewidth',2.5);

clear J_bc TP_bc
load('PF.mat');
for i = 1:N_UE
%     y = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded,3));
%     TP_bc(i) = y(end)/LTE_params.Tsubframe/1e6;
    TP_bc(i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
end
J_bc = sum(TP_bc)^2/(N_UE*sum(TP_bc.^2));
plot(J_bc,sum(TP_bc),'kd','Markersize',16,'linewidth',2.5);
% figure(2)
% plot(TP_bc(1),TP_bc(2),'kd','Markersize',14,'linewidth',2.5);

clear J_bc TP_bc
load('RF.mat');
LTE_params.scheduler.av_window = 500;
    a = [1,-(1-1/LTE_params.scheduler.av_window)];
    b = 1/LTE_params.scheduler.av_window;
for i = 1:N_UE
%     y = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded,3));
%     TP_bc(i) = y(end)/LTE_params.Tsubframe/1e6;
    TP_bc(i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
end
J_bc = sum(TP_bc)^2/(N_UE*sum(TP_bc.^2));
plot(J_bc,sum(TP_bc),'kv','Markersize',16,'linewidth',2.5);

clear J_bc TP_bc
load('RR.mat');
LTE_params.scheduler.av_window = 500;
    a = [1,-(1-1/LTE_params.scheduler.av_window)];
    b = 1/LTE_params.scheduler.av_window;
for i = 1:N_UE
%     y = filter(b,a,sum(simulation_results.UE_specific(i).throughput_coded,3));
%     TP_bc(i) = y(end)/LTE_params.Tsubframe/1e6;
    TP_bc(i) = mean((sum(simulation_results.UE_specific(i).throughput_coded,3))/LTE_params.Tsubframe/1e6);
end
J_bc = sum(TP_bc)^2/(N_UE*sum(TP_bc.^2));
plot(J_bc,sum(TP_bc),'k^','Markersize',16,'linewidth',2.5);

legend('SOCP','\alpha-fair','Max min','Best CQI','Prop. fair','Res. fair','Round robin');
