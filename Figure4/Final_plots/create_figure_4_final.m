%% Global settings
matlab_default_figure; %Changes the default figure appearance

% code_dir = '/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code';
% data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat';
% NSFR_fit_file = '/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data.mat';
% NSFR_sim_file = '/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data_predfr_all.mat';
% PLDS_sim_file = '/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_PLDS_data_predfr_all.mat';
% 
% cd(code_dir);

DATA_COLOR = [0 0 0];
NSFR_COLOR = [1 0 0];
PLDS_COLOR = [.4 .4 .4];

%% Panel A - Example data, plus NSFR and PLDS fits

%Make the figure;
load('Figure4/Final_plots/Panel_A');

neurons_to_show = neuron_stationarity_sort(5:-1:1);
fig1 = figure(11); clf
hold on;
for i1 = 1:length(neurons_to_show)
  plot(resps_orig(neurons_to_show(i1),:)'*20+(i1-1)*20, 'Color', DATA_COLOR);
  line([0,100],[(i1-1)*20,(i1-1)*20],'Color','k');
end

for i1 = 1:length(neurons_to_show)
  plot(fr_mean_pred_NSFR(neurons_to_show(i1),:)'*20+(i1-1)*20,'Color', NSFR_COLOR);
  scatter(params.ind_test, fr_mean_pred_NSFR(neurons_to_show(i1),params.ind_test)'*20+(i1-1)*20, 150,  NSFR_COLOR, 'x', 'LineWidth',2);
end
 

for i1 = 1:length(neurons_to_show)
  plot(fr_mean_pred_PLDS(neurons_to_show(i1),:)'*20+(i1-1)*20,'Color',PLDS_COLOR);
  scatter(params.ind_test, fr_mean_pred_PLDS(neurons_to_show(i1),params.ind_test)'*20+(i1-1)*20, 150,  PLDS_COLOR, 'x', 'LineWidth',2);
end

set(gca,'XTick',0:25:100);
ylim([0,20*length(neurons_to_show)]);
set(gca,'YTick',0:5:(20*length(neurons_to_show)-1));
set(gca,'YTickLabel',repmat(0:5:15,1,length(neurons_to_show)));
xlabel('Trial')
ylabel('Mean firing rate (Hz)');

set(fig1,'Units','centimeters');
set(fig1,'Position',[1,1,11,14]);

set(fig1,'PaperPositionMode','auto')
% print(fig1, '-depsc2', 'Figure4/Final_plots/Panel_A1.eps');
export_fig('Figure4/Final_plots/Panel_A1.eps', '-native');

neurons_to_show = neuron_stationarity_sort(58:62);
const_fr_step = 2;
fig2 = figure(12); clf
hold on;
for i1 = 1:length(neurons_to_show)
  plot(resps_orig(neurons_to_show(i1),:)'*20+(i1-1)*const_fr_step, 'Color', DATA_COLOR);
  line([0,100],[(i1-1)*const_fr_step,(i1-1)*const_fr_step],'Color','k');
end

for i1 = 1:length(neurons_to_show)
  plot(fr_mean_pred_NSFR(neurons_to_show(i1),:)'*20+(i1-1)*const_fr_step,'Color', NSFR_COLOR);
  scatter(params.ind_test, fr_mean_pred_NSFR(neurons_to_show(i1),params.ind_test)'*20+(i1-1)*const_fr_step, 100,  NSFR_COLOR, 'x', 'LineWidth',2);
end
 

for i1 = 1:length(neurons_to_show)
  plot(fr_mean_pred_PLDS(neurons_to_show(i1),:)'*20+(i1-1)*const_fr_step,'Color',PLDS_COLOR);
  scatter(params.ind_test, fr_mean_pred_PLDS(neurons_to_show(i1),params.ind_test)'*20+(i1-1)*const_fr_step, 100,  PLDS_COLOR, 'x', 'LineWidth',2);
end

set(gca,'XTick',0:25:100);
ylim([0,const_fr_step*length(neurons_to_show)]);
set(gca,'YTick',0:round(const_fr_step/2):(const_fr_step*length(neurons_to_show)-1));
set(gca,'YTickLabel',repmat(0:round(const_fr_step/2):(const_fr_step-1),1,length(neurons_to_show)));
xlabel('Trial')
ylabel('Mean firing rate (Hz)');

set(fig2,'Units','centimeters');
set(fig2,'Position',[1,1,11,14]);

set(fig2,'PaperPositionMode','auto')
% print(fig2, '-depsc2', 'Figure4/Final_plots/Panel_A2.eps');
export_fig('Figure4/Final_plots/Panel_A2.eps', '-native');

%% Panel B - Error values
load('Figure4/Final_plots/Panel_B.mat')
fig5 = figure(5); 
errorbar((1:8)-.05, mean(NSFR_error,2)/640*20, std(NSFR_error,0,2)/640*20, 'rx'); %/640 due to 10 left out trials and 64 neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
hold on; errorbar((1:8)+.05, mean(PLDS_error,2)/640*20, std(PLDS_error,0,2)/640*20, 'bx');
ylim([0.04 0.09])
xlabel('k')
set(gca,'XTick',1:8);
ylabel('RMSE'); %Population mean firing rate 4.5412 over everything;
legend('N-PLDS','PLDS')
set(gca,'YTick', 0.04:.01:0.09);
set(fig5,'Units','centimeters');
set(fig5,'Position',[1,1,25,13]);

set(fig5,'PaperPositionMode','auto')
% print(fig5, '-depsc2', 'Figure4/Final_plots/Panel_Btop.eps');
export_fig('Figure4/Final_plots/Panel_Btop.eps', '-native');

%Most non-stationary 5 neurons
neurons_to_include = neuron_stationarity_sort(1:5); %1 - most nonstationary, end - most stationary
NSFR_error = mean(NSFR_error_each(:,:,neurons_to_include),3);
PLDS_error = mean(PLDS_error_each(:,:,neurons_to_include),3);

fig6 = figure(6); 
errorbar((1:8)-.05, mean(NSFR_error,2)/(10*length(neurons_to_include))*20, std(NSFR_error,0,2)/(10*length(neurons_to_include))*20, 'rx'); %/(10*length(neurons_to_include)) due to 10 left out trials and certain number of chosen neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
hold on; errorbar((1:8)+.05, mean(PLDS_error,2)/(10*length(neurons_to_include))*20, std(PLDS_error,0,2)/(10*length(neurons_to_include))*20, 'bx');
xlabel('k')
set(gca,'XTick',1:8);
ylabel('RMSE'); %Population mean firing rate 4.5412 over everything;
legend('N-PLDS','PLDS')
set(fig6,'Units','centimeters');
set(fig6,'Position',[1,1,11,8]);

set(fig6,'PaperPositionMode','auto')
% print(fig6, '-depsc2', 'Figure4/Final_plots/Panel_Bbot1.eps');
export_fig('Figure4/Final_plots/Panel_Bbot1.eps', '-native');

%Most stationary 5 neurons
neurons_to_include = neuron_stationarity_sort(end-4:end); %1 - most nonstationary, end - most stationary
NSFR_error = mean(NSFR_error_each(:,:,neurons_to_include),3);
PLDS_error = mean(PLDS_error_each(:,:,neurons_to_include),3);

fig7 = figure(7); 
errorbar((1:8)-.05, mean(NSFR_error,2)/(10*length(neurons_to_include))*20, std(NSFR_error,0,2)/(10*length(neurons_to_include))*20, 'rx'); %/(10*length(neurons_to_include)) due to 10 left out trials and certain number of chosen neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
hold on; errorbar((1:8)+.05, mean(PLDS_error,2)/(10*length(neurons_to_include))*20, std(PLDS_error,0,2)/(10*length(neurons_to_include))*20, 'bx');
xlabel('k')
set(gca,'XTick',1:8);
ylabel('RMSE'); %Population mean firing rate 4.5412 over everything;
legend('N-PLDS','PLDS')
set(fig7,'Units','centimeters');
set(fig7,'Position',[1,1,11,8]);
set(fig7,'PaperPositionMode','auto')
% print(fig7, '-depsc2', 'Figure4/Final_plots/Panel_Bbot2.eps');
export_fig('Figure4/Final_plots/Panel_Bbot2.eps', '-native', fig7);

%% Panel C - Tau histogram
load('Figure4/Final_plots/Panel_C','tau2_end'); % k x runs tau^2 after last iteration
fig3 = figure(3);
tau2_end = tau2_end(:);
hist(sqrt(tau2_end),0:1:20); %Create histogram
xlim([0,20])
ylim([0,25]);
hold on;
mean_tau = mean(sqrt(tau2_end));
line([mean_tau mean_tau], [0,22],'Color','r','LineWidth',3); %Add line to the mean

xlabel('\tau (trials)');
ylabel('Counts');
% title('Tau is conserved over different latent dimensionalities and training sets');
set(fig3,'Units','centimeters');
set(fig3,'Position',[1,1,11,10]);
set(fig3,'PaperPositionMode','auto')
%print(fig3, '-depsc2', 'Figure4/Final_plots/Panel_C.eps');
export_fig('Figure4/Final_plots/Panel_C.eps', '-native', fig3);

%% Panel D - Total Covariance
load('Figure4/Final_plots/Panel_D');

%All XCov*_total are cell by timeshifts
%All XCov*_condi are trial by cell by timeshifts

%Compute normalized covariances
to_plot(1,:) = squeeze(mean(XCov_data_total,1))/max(mean(XCov_data_total,1));
to_plot(2,:) = squeeze(mean(XCov_NSFR_total,1))/max(mean(XCov_NSFR_total,1));
to_plot(3,:) = squeeze(mean(XCov_PLDS_total,1))/max(mean(XCov_PLDS_total,1));
to_plot(4,:) = squeeze(mean(mean(XCov_data_condi,1),2))/max(mean(mean(XCov_data_condi,1),2));
to_plot(5,:) = squeeze(mean(mean(XCov_NSFR_condi,1),2))/max(mean(mean(XCov_NSFR_condi,1),2));
to_plot(6,:) = squeeze(mean(mean(XCov_PLDS_condi,1),2))/max(mean(mean(XCov_PLDS_condi,1),2));


fig4 = figure(4); hold off;
%Total covariances
plot(-500:50:500,to_plot(1,:),'Color',DATA_COLOR);
hold on;
plot(-500:50:500,to_plot(2,:),'Color',NSFR_COLOR)
plot(-500:50:500,to_plot(3,:),'Color',PLDS_COLOR)

% %Conditional covariances
% plot(-500:50:500,to_plot(4,:),'Color',DATA_COLOR,'LineStyle','--'); 
% plot(-500:50:500,to_plot(5,:),'Color',NSFR_COLOR,'LineStyle','--')
% plot(-500:50:500,to_plot(6,:),'Color',PLDS_COLOR,'LineStyle','--')

h_d_leg = legend({'Data','N-PLDS (k=7)','PLDS (k=7)'});
h_d_xlab = xlabel('Time lag (ms)');
h_d_ylab = ylabel('Mean autocovariance (a.u.)');
set(h_d_leg, 'FontSize', 10)
set(fig4,'Units','centimeters');
set(fig4,'Position',[1,1,11,10]);
set(fig4,'PaperPositionMode','auto')
% print(fig4, '-depsc2', 'Figure4/Final_plots/Panel_D.eps');
export_fig('Figure4/Final_plots/Panel_D.eps', '-native');
