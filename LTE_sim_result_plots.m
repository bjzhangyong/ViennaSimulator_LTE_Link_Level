function LTE_sim_result_plots(simulation_results)
% Plots all the meaningful data from the results file
% Josep Colom Ikuno, jcolom@nt.tuwien.ac.at
% (c) 2009 by INTHFT
% www.nt.tuwien.ac.at

global LTE_params;

% Where to plot every figure
cell_BER_plot_figure        = 1;
cell_FER_plot_figure        = 2;
cell_throughput_plot_figure = 3;
user_BER_figure             = 4;
user_FER_figure             = 5;
user_throughput_figure      = 6;
user_BLER_figure            = 7;
MSE_freq_offset             = 9;

SNR_vector = simulation_results.SNR_vector;
minimum_for_semilogy = 10^-4;
confidence_intervals_color =  [0.8 0.8 0.8];

if LTE_params.plot_confidence_intervals
    %% Cell BER plot
   
    figure(cell_BER_plot_figure);
    semilogy(SNR_vector,simulation_results.cell_specific.BER_coded_overall,'b');
    hold on
    semilogy(SNR_vector,simulation_results.cell_specific.BER_uncoded_overall,'r');
    LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_BER_plot_figure,'ber', sum(simulation_results.cell_specific.biterrors_coded,3),sum(simulation_results.cell_specific.blocksize_coded,3));
    LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_BER_plot_figure,'ber', sum(simulation_results.cell_specific.biterrors_uncoded,3),sum(simulation_results.cell_specific.blocksize_uncoded,3));
    legend('cell coded BER','cell uncoded BER');
    xlabel('SNR [dB]');
    ylabel('BER');
    title('Cell BER');
    
    axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
    hold off
    grid on
    
    %% Cell FER/BLER plot
    figure(cell_FER_plot_figure);
    if isnan(simulation_results.cell_specific.FER_coded(:,:,2))
        semilogy(SNR_vector,mean(simulation_results.cell_specific.FER_coded(:,:,1),1),'b');
    else
        semilogy(SNR_vector,mean(mean(simulation_results.cell_specific.FER_coded,1),3),'b');
    end
    hold on
    if isnan(simulation_results.cell_specific.FER_uncoded(:,:,2))
        semilogy(SNR_vector,mean(simulation_results.cell_specific.FER_uncoded(:,:,1),1),'r');
    else
        semilogy(SNR_vector,mean(mean(simulation_results.cell_specific.FER_uncoded,1),3),'r');
    end
    
    if isnan(simulation_results.cell_specific.FER_coded(:,:,2))
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_FER_plot_figure,'mean', simulation_results.cell_specific.FER_coded(:,:,1));
    else
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_FER_plot_figure,'mean', mean(simulation_results.cell_specific.FER_coded,3));
    end
    if isnan(simulation_results.cell_specific.FER_uncoded(:,:,2))
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_FER_plot_figure,'mean', simulation_results.cell_specific.FER_uncoded(:,:,1));
    else
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_FER_plot_figure,'mean', mean(simulation_results.cell_specific.FER_uncoded,3));
    end
    legend('cell coded FER (BLER)','cell uncoded FER')
    xlabel('SNR [dB]')
    ylabel('FER')
    axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
    hold off
    grid on
    
    %% Cell throughtput plot
    figure(cell_throughput_plot_figure);
    
    % Plot total throughput (sum of all streams)
    cell_throughput_coded = sum(simulation_results.cell_specific.throughput_coded,3)/LTE_params.Tsubframe/1e6;
    plot(SNR_vector,mean(cell_throughput_coded,1),'.-b','Markersize',5);
    hold on
    cell_throughput_uncoded = sum(simulation_results.cell_specific.throughput_uncoded,3)/LTE_params.Tsubframe/1e6;
    plot(SNR_vector,mean(cell_throughput_uncoded,1),'.-r','Markersize',5);
    
    LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_throughput_plot_figure,'mean', cell_throughput_coded);
    LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,cell_throughput_plot_figure,'mean', cell_throughput_uncoded);
    
    legend('cell coded throughput','cell uncoded throughput','Location','best');
    xlabel('SNR [dB]');
    ylabel('Throughput [Mbit/s]');
    title('Cell throughput');
    hold off
    grid on
    
    %% User BER plot
    figure(user_BER_figure);
    BER_plot_matrix = zeros(length(SNR_vector),2*LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current_coded = sprintf('UE %d, coded',u_);
        current_uncoded = sprintf('UE %d, uncoded',u_);
        BER_plot_matrix(:,u_*2-1) = simulation_results.UE_specific(u_).BER_coded_overall';
        BER_plot_matrix(:,u_*2)   = simulation_results.UE_specific(u_).BER_uncoded_overall';
        legend_names_BER{u_*2-1} = current_coded;
        legend_names_BER{u_*2}   = current_uncoded;
    end
    
    semilogy(SNR_vector,BER_plot_matrix);
    hold on
    % which color??
    for u_ = 1:LTE_params.nUE
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_BER_figure,'ber', sum(simulation_results.UE_specific(u_).biterrors_coded,3),sum(simulation_results.UE_specific(u_).blocksize_coded,3));
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_BER_figure,'ber', sum(simulation_results.UE_specific(u_).biterrors_uncoded,3),sum(simulation_results.UE_specific(u_).blocksize_uncoded,3));
    end
    
    legend(legend_names_BER,'Location','best');
    xlabel('SNR [dB]');
    ylabel('BER');
    title('UE BER');
    axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
    hold off
    grid on
    
    %% User FER plot
    figure(user_FER_figure);
    if isnan(simulation_results.UE_specific(1).FER_coded(:,:,2))
        semilogy(SNR_vector,squeeze(mean(simulation_results.UE_specific(1).FER_coded(:,:,1),1)),'b');
    else
        semilogy(SNR_vector,mean(mean(simulation_results.UE_specific(1).FER_coded,1),3),'b');
    end
    hold on
    if isnan(simulation_results.UE_specific(1).FER_uncoded(:,:,2))
        semilogy(SNR_vector,squeeze(mean(simulation_results.UE_specific(1).FER_uncoded(:,:,1),1)),'r');
    else
        semilogy(SNR_vector,mean(mean(simulation_results.UE_specific(1).FER_uncoded,1),3),'r');
    end
    
    if isnan(simulation_results.cell_specific.FER_coded(:,:,2))
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_FER_figure,'mean', simulation_results.UE_specific(1).FER_coded(:,:,1));
    else
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_FER_figure,'mean', mean(simulation_results.UE_specific(1).FER_coded,3));
    end
    if isnan(simulation_results.cell_specific.FER_uncoded(:,:,2))
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_FER_figure,'mean', simulation_results.UE_specific(1).FER_uncoded(:,:,1));
    else
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_FER_figure,'mean', mean(simulation_results.UE_specific(1).FER_uncoded,3));
    end
    legend('first user coded FER','first user uncoded FER')
    xlabel('SNR [dB]')
    ylabel('FER')
    axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
    hold off
    grid on
    
    %% User throughtput plot
    figure(user_throughput_figure);
    throughput_plot_matrix = zeros(length(SNR_vector),2*LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current_coded = sprintf('UE %d, coded',u_);
        current_uncoded = sprintf('UE %d, uncoded',u_);
        throughput_plot_matrix(:,u_*2-1) = mean(sum(simulation_results.UE_specific(u_).throughput_coded,3))';
        throughput_plot_matrix(:,u_*2)   = mean(sum(simulation_results.UE_specific(u_).throughput_uncoded,3))';
        legend_names_throughput{u_*2-1} = current_coded;
        legend_names_throughput{u_*2}   = current_uncoded;
    end
    
    plot(SNR_vector,throughput_plot_matrix/LTE_params.Tsubframe/1e6);
    hold on
    for u_ = 1:LTE_params.nUE
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_throughput_figure,'mean', sum(simulation_results.UE_specific(u_).throughput_coded,3)/LTE_params.Tsubframe/1e6);
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_throughput_figure,'mean', sum(simulation_results.UE_specific(u_).throughput_uncoded,3)/LTE_params.Tsubframe/1e6);
    end
    
    legend(legend_names_throughput,'Location','best');
    xlabel('SNR [dB]');
    ylabel('Throughput [Mbit/s]');
    title('UE throughput');
    hold off
    grid on
    
    %% User BLER plot
    figure(user_BLER_figure);
    BLER_plot_matrix = zeros(length(SNR_vector),LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current = sprintf('UE %d',u_);
        BLER_plot_matrix(:,u_) = simulation_results.UE_specific(u_).BLER_overall';
        legend_names_BLER{u_} = current;
    end
    
    semilogy(SNR_vector,BLER_plot_matrix);
    hold on
    
    for u_ = 1:LTE_params.nUE
        LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_BLER_figure,'ber', sum(simulation_results.UE_specific(u_).FER_coded,3),sum(simulation_results.UE_specific(u_).used_codewords,3));
        %LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,user_BLER_figure,'ber', sum(simulation_results.UE_specific(u_).FER_uncoded,3),sum(simulation_results.UE_specific(u_).used_codewords,3));
    end
    legend(legend_names_BLER,'Location','best');
    
    plot_Y = simulation_results.UE_specific(1).BLER;
    plot_X = SNR_vector;
    
    grid on;
    xlabel('SNR [dB]');
    ylabel('BLER');
    title('BLER');
    axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
    legend(legend_names_BLER,'Location','best');
    
    %% User MSE carrier frequency offset
    figure(MSE_freq_offset)
    MSE_freq_offset_plot_matrix = zeros(length(SNR_vector),LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current = sprintf('UE %d',u_);
        MSE_freq_offset_plot_matrix(:,u_) = simulation_results.UE_specific(u_).MSE_freq_offset';
        legend_names_MSE_freq_offset{u_} = current;
    end

    semilogy(SNR_vector,MSE_freq_offset_plot_matrix);
    hold on
    
    if LTE_params.introduce_frequency_offset
        for u_ = 1:LTE_params.nUE
            LTE_plot_confidence_intervals(SNR_vector, confidence_intervals_color ,MSE_freq_offset,'mean', simulation_results.UE_specific(u_).freq_offset_est_error);
            %LTE_plot_confidence_intervals(SNR_vector, 'r' ,cell_BER_plot_figure,'mean', simulation_results.UE_specific(u_).freq_offset_est_error);
        end
    end
    
    legend(legend_names_MSE_freq_offset(u_),'Location','best');
    xlabel('SNR [dB]');
    ylabel('Mean Square Error');
    title('UE mean square error of the carrier frquency offset')
    hold off
    axis([min(SNR_vector) max(SNR_vector),minimum_for_semilogy 10^0])
    grid on
    grid on;
else
    %% Cell BER plot
    figure(cell_BER_plot_figure);
    semilogy(SNR_vector,simulation_results.cell_specific.BER_coded_overall,'b');
    hold on
    semilogy(SNR_vector,simulation_results.cell_specific.BER_uncoded_overall,'r');
    legend('cell coded BER','cell uncoded BER');
    xlabel('SNR [dB]');
    ylabel('BER');
    title('Cell BER');
    hold off
    grid on
    
    %% Cell FER/BLER plot
    % figure(cell_FER_plot_figure);
    % if isnan(simulation_results.cell_specific.FER_coded(:,:,2))
    %     semilogy(SNR_vector,mean(simulation_results.cell_specific.FER_coded(:,:,1),1),'b');
    % else
    %     semilogy(SNR_vector,mean(mean(simulation_results.cell_specific.FER_coded,1),3),'b');
    % end
    % hold on
    % if isnan(simulation_results.cell_specific.FER_uncoded(:,:,2))
    %     semilogy(SNR_vector,mean(simulation_results.cell_specific.FER_uncoded(:,:,1),1),'r');
    % else
    %     semilogy(SNR_vector,mean(mean(simulation_results.cell_specific.FER_uncoded,1),3),'r');
    % end
    % legend('cell coded FER (BLER)','cell uncoded FER')
    % xlabel('SNR [dB]')
    % ylabel('FER')
    % hold off
    % grid on
    
    %% Cell throughtput plot
    figure(cell_throughput_plot_figure);
    
    % Plot total throughput (sum of all streams)
    plot(SNR_vector,mean(sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6,'b'); 
    hold on 
    plot(SNR_vector,mean(sum(simulation_results.cell_specific.throughput_uncoded,3))/LTE_params.Tsubframe/1e6,'r');
    legend('cell coded throughput','cell uncoded throughput','Location','best');
    xlabel('SNR [dB]');
    ylabel('Throughput [Mbit/s]');
    title('Cell throughput');
    hold off
    grid on
    
    %% User BER plot
    figure(user_BER_figure);
    BER_plot_matrix = zeros(length(SNR_vector),2*LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current_coded = sprintf('UE %d, coded',u_);
        current_uncoded = sprintf('UE %d, uncoded',u_);
        BER_plot_matrix(:,u_*2-1) = simulation_results.UE_specific(u_).BER_coded_overall';
        BER_plot_matrix(:,u_*2)   = simulation_results.UE_specific(u_).BER_uncoded_overall';
        legend_names_BER{u_*2-1} = current_coded;
        legend_names_BER{u_*2}   = current_uncoded;
    end
    
    semilogy(SNR_vector,BER_plot_matrix);
    hold on
    legend(legend_names_BER,'Location','best');
    xlabel('SNR [dB]');
    ylabel('BER');
    title('UE BER');
    hold off
    grid on
    
    % %% User FER plot
    % figure(user_FER_figure);
    % if isnan(simulation_results.UE_specific(1).FER_coded(:,:,2))
    %     semilogy(SNR_vector,squeeze(mean(simulation_results.UE_specific(1).FER_coded(:,:,1),1)),'b');
    % else
    %     semilogy(SNR_vector,mean(mean(simulation_results.UE_specific(1).FER_coded,1),3),'b');
    % end
    % hold on
    % if isnan(simulation_results.UE_specific(1).FER_uncoded(:,:,2))
    %     semilogy(SNR_vector,squeeze(mean(simulation_results.UE_specific(1).FER_uncoded(:,:,1),1)),'r');
    % else
    %     semilogy(SNR_vector,mean(mean(simulation_results.UE_specific(1).FER_uncoded,1),3),'r');
    % end
    % legend('first user coded FER','first user uncoded FER')
    % xlabel('SNR [dB]')
    % ylabel('FER')
    % hold off
    % grid on
    
    %% User throughtput plot
    figure(user_throughput_figure);
    throughput_plot_matrix = zeros(length(SNR_vector),2*LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current_coded = sprintf('UE %d, coded',u_);
        current_uncoded = sprintf('UE %d, uncoded',u_);
        throughput_plot_matrix(:,u_*2-1) = mean(sum(simulation_results.UE_specific(u_).throughput_coded,3))';
        throughput_plot_matrix(:,u_*2)   = mean(sum(simulation_results.UE_specific(u_).throughput_uncoded,3))';
        legend_names_throughput{u_*2-1} = current_coded;
        legend_names_throughput{u_*2}   = current_uncoded;
    end
    
    plot(SNR_vector,throughput_plot_matrix/LTE_params.Tsubframe/1e6);
    hold on
    legend(legend_names_throughput,'Location','best');
    xlabel('SNR [dB]');
    ylabel('Throughput [Mbit/s]');
    title('UE throughput');
    hold off
    grid on
    
    %% User BLER plot
    figure(user_BLER_figure);
    BLER_plot_matrix = zeros(length(SNR_vector),LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current = sprintf('UE %d',u_);
        BLER_plot_matrix(:,u_) = simulation_results.UE_specific(u_).BLER_overall';
        legend_names_BLER{u_} = current;
    end
    
    semilogy(SNR_vector,BLER_plot_matrix);
    legend(legend_names_BLER,'Location','best');
    
    plot_Y = simulation_results.UE_specific(1).BLER;
    plot_X = SNR_vector;
    
    grid on;
    xlabel('SNR [dB]');
    ylabel('BLER');
    title('BLER');
    legend(legend_names_BLER,'Location','best');
    
    %% User MSE carrier frequency offset
    figure(MSE_freq_offset)
    MSE_freq_offset_plot_matrix = zeros(length(SNR_vector),LTE_params.nUE);
    
    for u_ = 1:LTE_params.nUE
        current = sprintf('UE %d',u_);
        MSE_freq_offset_plot_matrix(:,u_) = simulation_results.UE_specific(u_).MSE_freq_offset';
        legend_names_MSE_freq_offset{u_} = current;
    end
    
    semilogy(SNR_vector,MSE_freq_offset_plot_matrix);
    hold on
    legend(legend_names_MSE_freq_offset(u_),'Location','best');
    xlabel('SNR [dB]');
    ylabel('Mean Square Error');
    title('UE mean square error of the carrier frquency offset')
    hold off
    grid on
    grid on;
end