data_dir = '/nfs/data3/gergo/Mijung/Figure4/Output_data_final';
data_file = '/nfs/data3/gergo/Mijung/Figure4/Input_data/alexdata_session2_org0_bin50.mat';
code_dir = '/nfs/nhome/live/gbohner/Dropbox/Gatsby/Random/Mijung_NPLDS/Code';
output_data_file = 'final_NSFR_data.mat';
output_params_file = 'init_data.mat';


Model =  'NSFR';

fr_pred_error = zeros(10,10); % # of ks time # of runs
tau_hist = zeros(10,10,11); % k x run x [initial iter]
for k = 10
    for runnum = 2
      try
          [k runnum]
          output_data_dir = [data_dir filesep 'k' num2str(k) '_run' num2str(runnum)];
          saved_data_file = [output_data_dir filesep output_data_file];
          saved_params_file = [output_data_dir filesep output_params_file];
          fr_pred_error(k, runnum+1) = make_predictions( data_file, saved_data_file, saved_params_file, code_dir, k, Model,1);
      catch 
        printf('Error\n')
      end    
    end
%     figure(k);
%     show_p_value_contour(data_dir, output_data_file, k);
end




% sig_min = hp_pos(find(hp_val_mean==min(hp_val_mean)),1);
% 
% inds = find(hp_pos(:,1)==sig_min);
% 
% % errorbar(hp_pos(inds,2),hp_val_mean(inds),hp_val_std(inds));
% % plot(hp_pos(inds,2),hp_val(inds,:));
% % surf(tau_vals, sig_vals, log(log(hp_val_mean_2d)));
% % plot3_errorbars_surf(tau_vals, sig_vals, hp_val_mean_2d, hp_val_std_2d);
% 
% plot3_errorbars_surf(hp_pos(:,1),hp_pos(:,2),hp_val_mean,hp_val_std);

%  for k = 1:10
% figure(5); hold on; scatter(repmat(k,1,10), fr_pred_error(k,:),55,'x')
% end
% xlim([0,5])
% xlabel('k')
% title('RMSE in firing rates on left out trials (10%)')
% ylabel('RMSE')
% title('nPLDS RMSE in firing rates on left out trials (10%)')
