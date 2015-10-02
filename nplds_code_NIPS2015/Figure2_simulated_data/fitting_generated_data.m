%% to fit the data for Figure 2

clear all;
close all;
clc;

addpath ../standardEM/
addpath ../core_functions/

%%
% load data
load all_NSFR.mat

r = 100;

% put all the true params
params.A = A;
params.B = B;
params.C = C;
params.d = d;
params.h = h;

params.x0 = x0;
params.V0 = V0;
params.Q = eye(k);

newdataset = cell(r,1);

for i=1:r
    newdataset{i} = xyzinpn{i};
end

xyzinpn = newdataset;

%% start fitting here

Model = 'NSFR';
datastruct = VBEM_PLDSnonstationary(xyzinpn, r, params, Model); 

save fitting_NPLDS.mat

%% sanity check

z1_est = zeros(r, T);
z2_est = zeros(r, T);

z1 = zeros(r, T);
z2 = zeros(r, T);

corr_z_est = zeros(r,1);
corr_z = zeros(r, 1);

for trial_to_check = 1:r
    
    CC = datastruct.Mstep.C;
    hh = datastruct.Mstep.h(:, trial_to_check);
    covhh = datastruct.Mstep.covh(:,:,trial_to_check);
    mu = datastruct.Estep{trial_to_check}.mumarg;
    
    Cmud = zeros(p, T);
    for t=1:T
        Cmud(:,t) = CC*(mu(:,t)+hh) + datastruct.Mstep.d;
    end
    
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

subplot(311); plot(1:r, corr_z, 'o-', 1:r, corr_z_est, 'o-'); set(gca, 'ylim', [0 1]);
subplot(312); plot(1:T, sum(z1), 'k', 1:T, sum(z1_est), 'r')
subplot(313); plot(1:T, sum(z2), 'k', 1:T, sum(z2_est), 'r')

%%



