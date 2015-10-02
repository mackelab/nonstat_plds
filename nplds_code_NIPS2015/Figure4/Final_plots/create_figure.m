% %Best run -> k7 run 5 (index 6)
% PLDS_COLOR = 'b';
% 
% %% Panel A
data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat';
load(data_file, 'resps');
resps = double(resps);
resps = squeeze(mean(resps,2));
pop_mean = mean(mean(mean(resps)))*20; %Average spike / bin * 20 (due to 50 ms bins) = average firing rate;
resps_orig = resps;

%Do smoothing then least square to find non-stationarity
for i1 = 1:size(resps,1)
  resps(i1,:) = smooth(resps(i1,:),'moving',15);
end
resps_smooth = resps;
resps = bsxfun(@minus, resps, mean(resps,2));
resps = sqrt(mean(resps.^2,2));
[~, ix] = sort(resps,'descend');
% 
% % %Do linear regression to find non-stationarity
% % lincoeffs = zeros(size(resps,1),1);
% % for i1 = 1:size(resps,1)
% %   p = polyfit(1:100,resps(i1,:),1);
% %   lincoeffs(i1) = p(1);
% % end
% % [~,ix] = sort(abs(lincoeffs),'descend');
% 
% neurons_to_show = ix(1:3);
% fig1 = figure(1); clf
% hold on;
% for i1 = 1:length(neurons_to_show)
%   plot(resps_orig(neurons_to_show(i1),:)'*20+(i1-1)*20, 'black');
%   line([0,100],[(i1-1)*20,(i1-1)*20],'Color','k');
% end
%   
% 
% % %Compute predicted mean firing rates for training set;
% % resps_fit = zeros(size(resps_orig));
% % load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_NSFR_data.mat','datastruct', 'params');
% % for i1 = 1:length(params.ind_train)
% %   resps_fit(:,params.ind_train(i1)) = mean(exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d)),2);
% % end
% % hold on; plot(params.ind_train, resps_fit(neurons_to_show,params.ind_train)','r');
% 
% %Compute sampled mean firing rates for all data;
% load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data_predfr_all.mat', 'fr_mean_pred_all','params','hall', 'hpred');
% fr_mean_pred = median(fr_mean_pred_all,3);
% for i1 = 1:length(neurons_to_show)
%   plot(fr_mean_pred(neurons_to_show(i1),:)'*20+(i1-1)*20,'r');
% end
% 
% load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_PLDS_data_predfr_all.mat', 'fr_mean_pred_all','params','hall', 'hpred');
% fr_mean_pred = median(fr_mean_pred_all,3);
% for i1 = 1:length(neurons_to_show)
%   plot(fr_mean_pred(neurons_to_show(i1),:)'*20+(i1-1)*20,'Color',PLDS_COLOR);
% end
% 
% set(gca,'XTick',0:25:100);
% ylim([0,20*length(neurons_to_show)]);
% set(gca,'YTick',0:5:(20*length(neurons_to_show)-1));
% set(gca,'YTickLabel',repmat(0:5:15,1,length(neurons_to_show)));
% xlabel('Trial')
% ylabel('Firing rate (Hz)');
% 
% set(fig1,'Units','centimeters');
% set(fig1,'Position',[1,1,12,24]);
% 
% print('-depsc2', 'Figure4/Final_plots/Panel_A.eps');
% 
% %% Panel B
% %See h-predictions (they're on the line for h);
% load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data_predfr_all.mat', 'fr_mean_pred_all','params','hall', 'hpred');
% load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data.mat', 'datastruct');
% fig2 = figure(2);
% % plot(params.ind_train, exp(datastruct.Mstep{end}.C*datastruct.Mstep{end}.h)*20);
% % ylabel('Firing rate modulation');% AFter projecting through C matrix and exponentiating
% plot(params.ind_train, (datastruct.Mstep{end}.C*datastruct.Mstep{end}.h) + log(20)) % AFter projecting through C matrix
% ylabel('Log firing rate modulation');
% % hold on; scatter(reshape(repmat(params.ind_test,1,64),1,[]), reshape((datastruct.Mstep{end}.C*hpred)',1,[]));
% xlabel('Trial')
% set(fig2,'Units','centimeters');
% set(fig2,'Position',[1,1,20,10]);
% print('-depsc2', 'Figure4/Final_plots/Panel_B.eps');
% 
% 
%% Panel C
% load('Figure4/Final_plots/Panel_C.mat')
% fig3 = figure(3); %plot(tau_vals,tau_p_vals);
% load('hists.mat')
% tau_hist = tau_hist(1:9,:,end);
% tau_hist = tau_hist(:);
% hist(sqrt(tau_hist),0:1:20);
% xlim([0,20])
% 
% xlabel('tau (trials)');
% ylabel('Counts');
% title('Tau is conserved over different latent dimensionalities and training sets');
% set(fig3,'Units','centimeters');
% set(fig3,'Position',[1,1,20,10]);
% print('-depsc2', 'Figure4/Final_plots/Panel_C.eps');

%% Panel D - Sampled
%  max_lag = 10; %50 ms bins -> -500 ms to 500 ms
% %All the data
% %True_XCov
% load(data_file, 'resps');
% resps = double(resps);
% resps = reshape(resps,64,[]);
% XCov = zeros(size(resps,1), 2*max_lag+1);
% for i1 = 1:size(resps,1)
%     XCov(i1,:) = xcov(resps(i1,:),max_lag);
% end
% XCov_data_total = XCov;
% 
% %Trial Conditioned;
% %Data
% load(data_file, 'resps');
% resps = double(resps);
% XCov = zeros(size(resps,3), size(resps,1), 2*max_lag+1);
% for i0 = 1:size(resps(3))
%   for i1 = 1:size(resps,1)
%       XCov(i0, i1,:) = xcov(resps(i1,:,i0),max_lag);
%   end
% end
% XCov_data_condi = XCov;
% 
% % Compute the covariances for all sampled trials (given a certain parameter
% % setting), then average over the samples
% load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data_predfr_all.mat', 'fr_mean_pred_all','params','hall', 'hpred', 'all_data');
% % 
% %First concatanate trials from a single sample, compute the
% %cross_covariances, then average over all samples
% XCov = zeros(size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
% for j0 = 1:length(all_data{1})
%   fprintf('NSFR total %d\n', j0);
%   resps = [];
%   for j1 = 1:length(all_data)
%     resps = [resps all_data{j1}{j0}.y];
%   end
%   for i1 = 1:size(resps,1)
%     XCov(i1,:,j0) = xcov(resps(i1,:),max_lag);
%   end
% end
% XCov = median(XCov,3);
% XCov_NSFR_total = XCov;
% 
% XCov = zeros(length(all_data), size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
% for j0 = 1:length(all_data{1})
%   fprintf('NSFR condi %d\n', j0);
%   for j1 = 1:length(all_data)
%     resps = all_data{j1}{j0}.y;
%     for i1 = 1:size(resps,1)
%       XCov(j1,i1,:,j0) = xcov(resps(i1,:),max_lag);
%     end
%   end
% end
% XCov = median(XCov,4);
% XCov_NSFR_condi = XCov;
% 
% % Compute the covariances for all sampled trials (given a certain parameter
% % setting), then average over the samples
% load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_PLDS_data_predfr_all.mat', 'fr_mean_pred_all','params','hall', 'hpred', 'all_data');
% 
% %First concatanate trials from a single sample, compute the
% %cross_covariances, then average over all samples
% XCov = zeros(size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
% for j0 = 1:length(all_data{1})
%   fprintf('PLDS total %d\n', j0);
%   resps = [];
%   for j1 = 1:length(all_data)
%     resps = [resps all_data{j1}{j0}.y];
%   end
%   for i1 = 1:size(resps,1)
%     XCov(i1,:,j0) = xcov(resps(i1,:),max_lag);
%   end
% end
% XCov = median(XCov,3);
% XCov_PLDS_total = XCov;
% 
% XCov = zeros(length(all_data), size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
% for j0 = 1:length(all_data{1})
%   fprintf('PLDS condi %d\n', j0);
%   for j1 = 1:length(all_data)
%     resps = all_data{j1}{j0}.y;
%     for i1 = 1:size(resps,1)
%       XCov(j1,i1,:,j0) = xcov(resps(i1,:),max_lag);
%     end
%   end
% end
% XCov = median(XCov,4);
% XCov_PLDS_condi = XCov;
% 
% save('Figure4/Final_plots/Panel_D','XCov_*') 

% load('Figure4/Final_plots/Panel_D');
% 
% to_plot(1,:) = squeeze(mean(XCov_data_total,1))/max(mean(XCov_data_total,1));
% to_plot(2,:) = squeeze(mean(XCov_NSFR_total,1))/max(mean(XCov_NSFR_total,1));
% to_plot(3,:) = squeeze(mean(XCov_PLDS_total,1))/max(mean(XCov_PLDS_total,1));
% % to_plot(4,:) = squeeze(mean(mean(XCov_data_condi,1),2))/max(mean(mean(XCov_data_condi,1),2));
% % to_plot(5,:) = squeeze(mean(mean(XCov_NSFR_condi,1),2))/max(mean(mean(XCov_NSFR_condi,1),2));
% % to_plot(6,:) = squeeze(mean(mean(XCov_PLDS_condi,1),2))/max(mean(mean(XCov_PLDS_condi,1),2));
% 
% % to_plot(1,:) = squeeze(mean(XCov_data_total,1));
% % to_plot(2,:) = squeeze(mean(XCov_NSFR_total,1));
% % to_plot(3,:) = squeeze(mean(XCov_PLDS_total,1));
% % to_plot(4,:) = squeeze(mean(mean(XCov_data_condi,1),2));
% % to_plot(5,:) = squeeze(mean(mean(XCov_NSFR_condi,1),2));
% % to_plot(6,:) = squeeze(mean(mean(XCov_PLDS_condi,1),2));
% 
% fig4 = figure(4); hold off;
% plot(-500:50:500,to_plot(1,:),'k');
% hold on;
% plot(-500:50:500,to_plot(2,:),'r')
% plot(-500:50:500,to_plot(3,:),'b')
% 
% % plot(-500:50:500,to_plot(4,:),'k--'); 
% % hold on;
% % plot(-500:50:500,to_plot(5,:),'r--')
% % plot(-500:50:500,to_plot(6,:),'b--')
% 
% legend({'Data','NSFR','PLDS'})
% xlabel('Time lag (ms)')
% ylabel('Normalized autocovariance')
% set(fig4,'Units','centimeters');
% set(fig4,'Position',[1,1,15,10]);
% print('-depsc2', 'Figure4/Final_plots/Panel_D_total.eps');

%% Panel D - Sampled neuron pairs
 max_lag = 10; %50 ms bins -> -500 ms to 500 ms
%All the data
%True_XCov
load(data_file, 'resps');
resps = double(resps);
resps = reshape(resps,64,[]);
XCov = zeros(size(resps,1), size(resps,1), 2*max_lag+1);
for i1 = 1:size(resps,1)
  for i2 = 1:size(resps,1)
    XCov(i1,i2,:) = xcov(resps(i1,:),resps(i2,:),max_lag);
  end
end
XCov_data_total = XCov;

%Trial Conditioned;
%Data
load(data_file, 'resps');
resps = double(resps);
XCov = zeros(size(resps,3), size(resps,1),size(resps,1), 2*max_lag+1);
for i0 = 1:size(resps(3))
  for i1 = 1:size(resps,1)
    for i2 = 1:size(resps,1)
      XCov(i0, i1,i2,:) = xcov(resps(i1,:),resps(i2,:),max_lag);
    end
  end
end
XCov_data_condi = XCov;

% Compute the covariances for all sampled trials (given a certain parameter
% setting), then average over the samples
load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_NSFR_data_predfr_all.mat', 'fr_mean_pred_all','params','hall', 'hpred', 'all_data');
% 
%First concatanate trials from a single sample, compute the
%cross_covariances, then average over all samples
XCov = zeros(size(resps,1), size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
for j0 = 1:length(all_data{1})
  fprintf('NSFR total %d\n', j0);
  resps = [];
  for j1 = 1:length(all_data)
    resps = [resps all_data{j1}{j0}.y];
  end
  for i1 = 1:size(resps,1)
    for i2 = 1:size(resps,1)
      XCov(i1,i2,:,j0) = xcov(resps(i1,:),resps(i2,:),max_lag);
    end
  end
end
XCov = median(XCov,4);
XCov_NSFR_total = XCov;

XCov = zeros(length(all_data), size(resps,1), size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
for j0 = 1:length(all_data{1})
  fprintf('NSFR condi %d\n', j0);
  for j1 = 1:length(all_data)
    resps = all_data{j1}{j0}.y;
    for i1 = 1:size(resps,1)
      for i2 = 1:size(resps,1)
        XCov(j1,i1, i2,:,j0) = xcov(resps(i1,:),resps(i2,:),max_lag);
      end
    end
  end
end
XCov = median(XCov,5);
XCov_NSFR_condi = XCov;

% Compute the covariances for all sampled trials (given a certain parameter
% setting), then average over the samples
load('/nfs/data3/gergo/Mijung/Figure4/Example_run/k7_run5/final_PLDS_data_predfr_all.mat', 'fr_mean_pred_all','params','hall', 'hpred', 'all_data');

%First concatanate trials from a single sample, compute the
%cross_covariances, then average over all samples
XCov = zeros(size(resps,1), size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
for j0 = 1:length(all_data{1})
  fprintf('PLDS total %d\n', j0);
  resps = [];
  for j1 = 1:length(all_data)
    resps = [resps all_data{j1}{j0}.y];
  end
  for i1 = 1:size(resps,1)
    for i2 = 1:size(resps,1)
      XCov(i1,i2,:,j0) = xcov(resps(i1,:),resps(i2,:),max_lag);
    end
  end
end
XCov = median(XCov,4);
XCov_PLDS_total = XCov;

XCov = zeros(length(all_data), size(resps,1), size(resps,1), 2*max_lag+1, length(all_data{1})); %neuron by lag by samples
for j0 = 1:length(all_data{1})
  fprintf('PLDS condi %d\n', j0);
  for j1 = 1:length(all_data)
    resps = all_data{j1}{j0}.y;
    for i1 = 1:size(resps,1)
      for i2 = 1:size(resps,1)
        XCov(j1,i1, i2,:,j0) = xcov(resps(i1,:),resps(i2,:),max_lag);
      end
    end
  end
end
XCov = median(XCov,5);
XCov_PLDS_condi = XCov;

save('Figure4/Final_plots/Panel_D_neur_pairs','XCov_*') 

load('Figure4/Final_plots/Panel_D_neur_pairs');

 to_plot(1,:) = squeeze(mean(mean(XCov_data_total,1),2))/max(mean(mean(XCov_data_total,1),2));
to_plot(2,:) = squeeze(mean(mean(XCov_NSFR_total,1),2))/max(mean(mean(XCov_NSFR_total,1),2));
to_plot(3,:) = squeeze(mean(mean(XCov_PLDS_total,1),2))/max(mean(mean(XCov_PLDS_total,1),2));
to_plot(4,:) = squeeze(mean(mean(mean(XCov_data_condi,1),2),3))/max(mean(mean(mean(XCov_data_condi,1),2),3));
to_plot(5,:) = squeeze(mean(mean(mean(XCov_NSFR_condi,1),2),3))/max(mean(mean(mean(XCov_NSFR_condi,1),2),3));
to_plot(6,:) = squeeze(mean(mean(mean(XCov_PLDS_condi,1),2),3))/max(mean(mean(mean(XCov_PLDS_condi,1),2),3));

fig4 = figure(4); hold off;
plot(-500:50:500,to_plot(1,:),'k');
hold on;
plot(-500:50:500,to_plot(2,:),'r')
plot(-500:50:500,to_plot(3,:),'b')

plot(-500:50:500,to_plot(4,:),'k--'); 
hold on;
plot(-500:50:500,to_plot(5,:),'r--')
plot(-500:50:500,to_plot(6,:),'b--')

legend({'Data','NSFR','PLDS'})
xlabel('Time lag (ms)')
ylabel('Normalized autocovariance')
set(fig4,'Units','centimeters');
set(fig4,'Position',[1,1,15,10]);
print('-depsc2', 'Figure4/Final_plots/Panel_D_neur_pairs.eps');

%% Panel D
% max_lag = 10; %50 ms bins -> -500 ms to 500 ms
% %All the data
% %True_XCov
% load(data_file, 'resps');
% resps = double(resps);
% resps = reshape(resps,64,[]);
% XCov = zeros(size(resps,1), size(resps,1), 2*max_lag+1);
% for i1 = 1:size(resps,1)
%   for i2 = 1:size(resps,1)
%     XCov(i1,i2,:) = xcov(resps(i1,:),resps(i2,:),max_lag);
%   end
% end
% XCov_data_total = XCov;
% 
% %NSFR_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_NSFR_data.mat','datastruct');
% resps = [];
% for i1 = 1:length(params.ind_train)
%   resps = [resps, exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d))];
% end
% XCov = zeros(size(resps,1), size(resps,1), 2*max_lag+1);
% for i1 = 1:size(resps,1)
%   for i2 = 1:size(resps,1)
%     XCov(i1,i2,:) = xcov(resps(i1,:),resps(i2,:),max_lag);
%   end
% end
% XCov_NSFR_total = XCov;
% 
% %PLDS_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_PLDS_data.mat','datastruct');
% resps = [];
% for i1 = 1:length(params.ind_train)
%   resps = [resps, exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d))];
% end
% XCov = zeros(size(resps,1), size(resps,1), 2*max_lag+1);
% for i1 = 1:size(resps,1)
%   for i2 = 1:size(resps,1)
%     XCov(i1,i2,:) = xcov(resps(i1,:),resps(i2,:),max_lag);
%   end
% end
% XCov_PLDS_total = XCov;
% % 
% % 
% % 
% % 
% %Trial Conditioned;
% %Data
% load(data_file, 'resps');
% resps = double(resps);
% XCov = zeros(size(resps,3), size(resps,1), size(resps,1), 2*max_lag+1);
% for i0 = 1:size(resps(3))
%   for i1 = 1:size(resps,1)
%     for i2 = 1:size(resps,1)
%       XCov(i0, i1,i2,:) = xcov(resps(i1,:,i0),resps(i2,:,i0),max_lag);
%     end
%   end
% end
% XCov_data_condi = XCov;
% 
% %NSFR_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_NSFR_data.mat','datastruct');
% resps = zeros(64,80,90);
% for i1 = 1:length(params.ind_train)
%   resps(:,:,i1) = exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d));
% end
% XCov = zeros(size(resps,3), size(resps,1), size(resps,1), 2*max_lag+1);
% for i0 = 1:size(resps(3))
%   for i1 = 1:size(resps,1)
%     for i2 = 1:size(resps,1)
%       XCov(i0, i1,i2,:) = xcov(resps(i1,:,i0),resps(i2,:,i0),max_lag);
%     end
%   end
% end
% XCov_NSFR_condi = XCov;
% 
% %PLDS_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_PLDS_data.mat','datastruct');
% resps = zeros(64,80,90);
% for i1 = 1:length(params.ind_train)
%   resps(:,:,i1) = exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d));
% end
% XCov = zeros(size(resps,3), size(resps,1), size(resps,1), 2*max_lag+1);
% for i0 = 1:size(resps(3))
%   for i1 = 1:size(resps,1)
%     for i2 = 1:size(resps,1)
%       XCov(i0, i1,i2,:) = xcov(resps(i1,:,i0),resps(i2,:,i0),max_lag);
%     end
%   end
% end
% XCov_PLDS_condi = XCov;
% % 
% save('Figure4/Final_plots/Panel_D','XCov_*') 

% load('Figure4/Final_plots/Panel_D');
% 
% to_plot(1,:) = squeeze(mean(mean(XCov_data_total,1),2))/max(mean(mean(XCov_data_total,1),2));
% to_plot(2,:) = squeeze(mean(mean(XCov_NSFR_total,1),2))/max(mean(mean(XCov_NSFR_total,1),2));
% to_plot(3,:) = squeeze(mean(mean(XCov_PLDS_total,1),2))/max(mean(mean(XCov_PLDS_total,1),2));
% to_plot(4,:) = squeeze(mean(mean(mean(XCov_data_condi,1),2),3))/max(mean(mean(mean(XCov_data_condi,1),2),3));
% to_plot(5,:) = squeeze(mean(mean(mean(XCov_NSFR_condi,1),2),3))/max(mean(mean(mean(XCov_NSFR_condi,1),2),3));
% to_plot(6,:) = squeeze(mean(mean(mean(XCov_PLDS_condi,1),2),3))/max(mean(mean(mean(XCov_PLDS_condi,1),2),3));
% 
% figure(4); hold off;
% plot(-500:50:500,to_plot(1,:),'k');
% hold on;
% plot(-500:50:500,to_plot(2,:),'r')
% plot(-500:50:500,to_plot(3,:),'b')
% 
% plot(-500:50:500,to_plot(4,:),'k--'); 
% hold on;
% plot(-500:50:500,to_plot(5,:),'r--')
% plot(-500:50:500,to_plot(6,:),'b--')

%% Panel D _ new
% max_lag = 10; %50 ms bins -> -500 ms to 500 ms
% %All the data
% %True_XCov
% load(data_file, 'resps');
% resps = double(resps);
% resps = reshape(resps,64,[]);
% XCov = zeros(size(resps,1), 2*max_lag+1);
% for i1 = 1:size(resps,1)
%     XCov(i1,:) = xcov(resps(i1,:),max_lag);
% end
% XCov_data_total = XCov;
% 
% %NSFR_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_NSFR_data.mat','datastruct');
% resps = [];
% for i1 = 1:length(params.ind_train)
%   resps = [resps, exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d))];
% end
% XCov = zeros(size(resps,1), 2*max_lag+1);
% for i1 = 1:size(resps,1)
%     XCov(i1,:) = xcov(resps(i1,:),max_lag);
% end
% XCov_NSFR_total = XCov;
% 
% %PLDS_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_PLDS_data.mat','datastruct');
% resps = [];
% for i1 = 1:length(params.ind_train)
%   resps = [resps, exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d))];
% end
% XCov = zeros(size(resps,1), 2*max_lag+1);
% for i1 = 1:size(resps,1)
%     XCov(i1,:) = xcov(resps(i1,:),max_lag);
% end
% XCov_PLDS_total = XCov;
% % 
% % 
% % 
% % 
% %Trial Conditioned;
% %Data
% load(data_file, 'resps');
% resps = double(resps);
% XCov = zeros(size(resps,3), size(resps,1), 2*max_lag+1);
% for i0 = 1:size(resps(3))
%   for i1 = 1:size(resps,1)
%       XCov(i0, i1,:) = xcov(resps(i1,:,i0),max_lag);
%   end
% end
% XCov_data_condi = XCov;
% 
% %NSFR_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_NSFR_data.mat','datastruct');
% resps = zeros(64,80,90);
% for i1 = 1:length(params.ind_train)
%   resps(:,:,i1) = exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d));
% end
% XCov = zeros(size(resps,3), size(resps,1), 2*max_lag+1);
% for i0 = 1:size(resps(3))
%   for i1 = 1:size(resps,1)
%       XCov(i0, i1,:) = xcov(resps(i1,:,i0),max_lag);
%   end
% end
% 
% XCov_NSFR_condi = XCov;
% 
% %PLDS_XCov
% load('/nfs/data3/gergo/Mijung/Figure4/Output_data_final/k7_run5/final_PLDS_data.mat','datastruct');
% resps = zeros(64,80,90);
% for i1 = 1:length(params.ind_train)
%   resps(:,:,i1) = exp(bsxfun(@plus,datastruct.Mstep{end}.C*(bsxfun(@plus,datastruct.Estep{end}{1}.mumarg,datastruct.Mstep{end}.h(:,i1))), datastruct.Mstep{end}.d));
% end
% XCov = zeros(size(resps,3), size(resps,1), 2*max_lag+1);
% for i0 = 1:size(resps(3))
%   for i1 = 1:size(resps,1)
%       XCov(i0, i1,:) = xcov(resps(i1,:,i0),max_lag);
%   end
% end
% 
% XCov_PLDS_condi = XCov;
% 
% save('Figure4/Final_plots/Panel_D_auto','XCov_*') 
% 
% load('Figure4/Final_plots/Panel_D_auto');
% 
% to_plot(1,:) = squeeze(mean(XCov_data_total,1))/max(mean(XCov_data_total,1));
% to_plot(2,:) = squeeze(mean(XCov_NSFR_total,1))/max(mean(XCov_NSFR_total,1));
% to_plot(3,:) = squeeze(mean(XCov_PLDS_total,1))/max(mean(XCov_PLDS_total,1));
% to_plot(4,:) = squeeze(mean(mean(XCov_data_condi,1),2))/max(mean(mean(XCov_data_condi,1),2));
% to_plot(5,:) = squeeze(mean(mean(XCov_NSFR_condi,1),2))/max(mean(mean(XCov_NSFR_condi,1),2));
% to_plot(6,:) = squeeze(mean(mean(XCov_PLDS_condi,1),2))/max(mean(mean(XCov_PLDS_condi,1),2));
% 
% figure(4); hold off;
% plot(-500:50:500,to_plot(1,:),'k');
% hold on;
% plot(-500:50:500,to_plot(2,:),'r')
% plot(-500:50:500,to_plot(3,:),'b')
% 
% plot(-500:50:500,to_plot(4,:),'k--'); 
% hold on;
% plot(-500:50:500,to_plot(5,:),'r--')
% plot(-500:50:500,to_plot(6,:),'b--')

%% Panel E
load('Figure4/Final_plots/Panel_E.mat')
fig5 = figure(5); 
errorbar((1:9)-.1, mean(NSFR_error,2)/640*20, std(NSFR_error,0,2)/640*20, 'rx'); %/640 due to 10 left out trials and 64 neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
hold on; errorbar((1:9)+.1, mean(PLDS_error,2)/640*20, std(PLDS_error,0,2)/640*20, 'bx');
ylim([0.002 0.005]*20)
xlabel('k')
title('RMSE in firing rates on left out trials (10%) after 30 iter')
%ylabel('RMSE')
ylabel('RMSE in predicted single neuron firing rate'); %Population mean 4.5412 (pop_mean);
legend('NSFR','PLDS')
set(fig5,'Units','centimeters');
set(fig5,'Position',[1,1,20,10]);
print('-depsc2', 'Figure4/Final_plots/Panel_E.eps');
% 
%% Panel F
% Compute prediction error for just the X most non-stationary neurons 
bad_neurs = ix(:); %Sorted from most non-stationary to most stationary

data_dir = '/nfs/data3/gergo/Mijung/Figure4/Output_data_final';
data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat';
code_dir = '/nfs/nhome/live/gbohner/Git_main/nonstat_plds_code/Gergo';
output_data_file_temp = 'final_MODEL_data.mat';
% output_data_file_temp = 'datastruct_MODEL_iter_10.mat';
output_params_file = 'init_data.mat';
% 
% clf;
% Model =  'PLDS';
Models = {'NSFR','PLDS'};

for mod = 1:2
Model = Models{mod};
output_data_file = strrep(output_data_file_temp, 'MODEL', Model);
fr_pred_error_each = zeros(9,10, 64); % # of ks by # of runs by neurons
for k = 1:9%
    for runnum = 0:9
       output_data_dir = [data_dir filesep 'k' num2str(k) '_run' num2str(runnum)];
        saved_data_file = [output_data_dir filesep output_data_file];
        saved_params_file = [output_data_dir filesep output_params_file];
        load([saved_data_file(1:end-4) '_predfr.mat'], 'Mstepresults', 'hpred', 'params','fr_mean_pred_all', 'fr_mean_true');
        fr_mean_pred = median(fr_mean_pred_all,3);
        fr_pred_error_each(k, runnum+1,:) = sqrt(sum((fr_mean_pred - fr_mean_true).^2,2));
    end
end


if mod == 1
NSFR_error_each = fr_pred_error_each;
else
PLDS_error_each = fr_pred_error_each;
end
% 
end
% 
% figure(6); hold off;
% errorbar((1:9)-.1, mean(NSFR_error,2)/(10*length(bad_neurs))*20, std(NSFR_error,0,2)/(10*length(bad_neurs))*20, 'rx'); %/640 due to 10 left out trials and 64 neurons/trials to get error / single bin, *20 to get back to Hz from 50 ms bins.
% hold on; errorbar((1:9)+.1, mean(PLDS_error,2)/(10*length(bad_neurs))*20, std(PLDS_error,0,2)/(10*length(bad_neurs))*20, 'bx');
% % ylim([0.000 0.02]*20)
% xlabel('k')
% title('RMSE in firing rates on left out trials (10%) after 30 iter')
% %ylabel('RMSE')
% ylabel('RMSE in predicted single neuron firing rate'); %Population mean 4.5412 (pop_mean);
% legend('NSFR','PLDS')