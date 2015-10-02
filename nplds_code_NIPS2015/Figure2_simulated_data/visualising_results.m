
clear all;
clc;
close all;

load fitting_NPLDS.mat;

z1_est = zeros(r, T);
z2_est = zeros(r, T);

z1 = zeros(r, T);
z2 = zeros(r, T);

corr_z_est = zeros(r,1);
corr_z = zeros(r, 1);

zestmat = zeros(p, T, r); 

for trial_to_check = 1:r
    
    CC = datastruct.Mstep.C;
    hh = datastruct.Mstep.h(:, trial_to_check);
    covhh = datastruct.Mstep.covh(:,:,trial_to_check);
    mu = datastruct.Estep{trial_to_check}.mumarg;
    
    Cmud = zeros(p, T);
    for t=1:T
        Cmud(:,t) = CC*(mu(:,t)+hh) + datastruct.Mstep.d;
    end
    
    zestmat(:,:,trial_to_check) = Cmud;
    
    z1_est(trial_to_check,:) = sum(Cmud(1:p/2,:));
    z2_est(trial_to_check,:) = sum(Cmud(p/2+1:p,:));
    
    corr_z_est(trial_to_check) = corr(z1_est(trial_to_check,:)', z2_est(trial_to_check,:)');
    
    invsig = datastruct.Estep{trial_to_check}.inv_sigmarg;
    covz_errbar = zeros(p, p, T);
    errorbar = zeros(p, T);
    for t=1:T
        covz_errbar(:,:,t) = CC*(inv(invsig(:,:,t))+covhh)*CC';
        errorbar(:,t) =  diag(covz_errbar(:,:,t));
    end
    
    
    % plot(1:T, k1k2mat(2,:,1)', 'k', 1:T, k1k2mat_est(2,:,1), 'r', 1:T, k1k2mat_est(2,:,1)-1.64*sqrt(sum(errorbar(1+p/2:p,:))), 'r--', 1:T, k1k2mat_est(2,:,1)+1.64*sqrt(sum(errorbar(1+p/2:p,:))), 'r--');
    % set(gca, 'ylim', [-45 -35])
    
    z1(trial_to_check,:) = sum(xyzinpn{trial_to_check}.z(1:p/2,:,1));
    z2(trial_to_check,:) = sum(xyzinpn{trial_to_check}.z(p/2+1:p,:,1));
    
    corr_z(trial_to_check) = corr(z1(trial_to_check,:)', z2(trial_to_check,:)'); 

    
end

%%
    

figure
subplot(211);
plot(1:r, mean(z1,2)/(p/2),'r', 1:r, mean(z2,2)/(p/2), 'b', ...
    1:r, mean(z1_est,2)/(p/2), 'r--', 1:r, mean(z2_est,2)/(p/2), 'b--')

%%
% subplot(311); plot(1:r, corr_z, 'o-', 1:r, corr_z_est, 'o-'); 
% set(gca, 'ylim', [0 1]);
% subplot(211); plot(1:T, sum(z1), 'k', 1:T, sum(z1_est), 'r')
% subplot(212); plot(1:T, sum(z2), 'k', 1:T, sum(z2_est), 'r')
% 
% alpha = 2.5;
% % alpha =10;
% subplot(211); plot(1:T, sum(z1), 'k', 1:T, sum(z1_est), 'r', 1:T, sum(z1_est)-alpha*sum(sqrt(errorbar(1:p/2,:))), 'r--', 1:T, sum(z1_est)+alpha*sum(sqrt(errorbar(1:p/2,:))), 'r--')
% subplot(212); plot(1:T, sum(z2), 'k', 1:T, sum(z2_est), 'r', 1:T, sum(z2_est)-alpha*sum(sqrt(errorbar(1+p/2:p,:))), 'r--', 1:T, sum(z2_est)+alpha*sum(sqrt(errorbar(1+p/2:p,:))), 'r--')

%%

numtimebins = T;
autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
%    autocorr(whichcell, :) = xcov(yy(whichcell,:), numtimebins/2, 'unbiased');
     autocorr(whichcell, :) = xcov(zz(whichcell,:), numtimebins/2, 'unbiased');
end

avgautocorr_acrcells = mean(autocorr);

% figure; plot(autocorr'); title('total auto-cov of each cell');

autocorr_condi = zeros(p, numtimebins+1, r);

for whichcell = 1:p
    
    for whichtrial = 1:r
        autocorr_condi(whichcell, :, whichtrial) = xcov(z(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
%        autocorr_condi(whichcell, :, whichtrial) = xcov(y(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
    end
    
end

autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

figure
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r')
legend('total covariance', 'conditional covariance');
%% 


zzz = [];
for i=1:r
    zzz = [zz zestmat(:,:,i)];
end

% numtimebins = T/4;
autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
%    autocorr(whichcell, :) = xcov(yy(whichcell,:), numtimebins/2, 'unbiased');
     autocorr(whichcell, :) = xcov(zzz(whichcell,:), numtimebins/2, 'unbiased');
end

avgautocorr_acrcells = mean(autocorr);

% figure; plot(autocorr'); title('total auto-cov of each cell');

autocorr_condi = zeros(p, numtimebins+1, r);

for whichcell = 1:p
    
    for whichtrial = 1:r
        autocorr_condi(whichcell, :, whichtrial) = xcov(zestmat(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
%        autocorr_condi(whichcell, :, whichtrial) = xcov(y(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
    end
    
end

autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

hold on; 
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k--', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r--')

%%

ztrue_acrsstrial = squeeze(mean(z, 2));
