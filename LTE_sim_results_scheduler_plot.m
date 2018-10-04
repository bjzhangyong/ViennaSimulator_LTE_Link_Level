function LTE_sim_results_scheduler_plot(BS,BS_output,TTI_num)
% Plots a plot of the scheduling. Including reference and sync signals.
% Author: Michal Simko, msimko@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

global LTE_params;

%% Scheduler plot
scheduler_figure = figure(8);
scheduler_freq = zeros(LTE_params.Ntot,2);
for u_id = 1:LTE_params.nUE
    scheduler_freq(logical(kron(BS_output.UE_signaling(u_id).MCS_and_scheduling.UE_mapping,ones(LTE_params.Nsc,1)))) = u_id;
end
scheduler_matrix = [repmat(scheduler_freq(:,1),1,LTE_params.Ns),repmat(scheduler_freq(:,2),1,LTE_params.Ns)];
subframe_corr = 1;
RefSym          = LTE_params.Reference_Signal(1,subframe_corr).RefSym;
RefMapping      = LTE_params.Reference_Signal(1,subframe_corr).RefMapping;
total_no_refsym = LTE_params.Reference_Signal(1,subframe_corr).total_no_refsym;
NoData_indices  = LTE_params.Reference_Signal(1,subframe_corr).NoData_indices;
for ii = 1:BS.nTX
    scheduler_matrix(squeeze(RefMapping(:,:,ii))) = ii+LTE_params.nUE;
end
scheduler_YTick = 0.5:12:LTE_params.Nrb*12+0.5;
scheduler_XTick = 0.5:2*LTE_params.Ns+0.5;
for i_=1:length(scheduler_YTick)
    scheduler_YTickLabel{i_} = num2str(scheduler_YTick(i_)-0.5);
end
for i_=1:length(scheduler_XTick)
    scheduler_XTickLabel{i_} = num2str(scheduler_XTick(i_)-0.5);
end

scheduler_axes = axes('Parent',scheduler_figure,...
    'YTick',scheduler_YTick,...
    'YTickLabel',scheduler_YTickLabel,...
    'YDir','reverse',...
    'Layer','top',...
    'XTick',scheduler_XTick,...
    'XTickLabel',scheduler_XTickLabel,...
    'Layer','top');
box('on');
grid('on');
hold('all');
imagesc(scheduler_matrix,'Parent',scheduler_axes,'CDataMapping','scaled');
xlim([0.5 2*LTE_params.Ns+0.5]);
ylim([0.5 LTE_params.Nrb*12+0.5]);
xlabel('time index (symbol)');
ylabel('frequency index (subcarrier number)');
%legend('user1','user2','pilot1','pilot2','pilot3','pilots4');  % Legend does not work here
title(sprintf('scheduler plot (TTI %3.0f)',TTI_num));
