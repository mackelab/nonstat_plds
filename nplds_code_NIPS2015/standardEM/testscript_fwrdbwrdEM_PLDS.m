%% To test EM with forward/backward algorithm 
% May 04, 2014
% written by Mijung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

% this is for hinton diagram
set(0, 'DefaultFigureMenuBar', 'figure');
colordef white;

addpath ../mattBealsCode_v3_4_1/

%% define essential quantities

% k << pinp and k<p
k = 2; % dimensionality of hidden states
pinp = 5; % input dimension (d, in my writeup)
p = 10; % output dimension (e.g., # of neurons)
T = 100; % time 1:T (for training)

%% generate params

trueparams = generate_params(k, pinp, p);

%% generate inputs (if pinp!=0)

if pinp>0
    trueparams.inpn = generate_inputs(T, pinp);
else
    trueparams.inpn = [];
end

%% generate latent variables/observations
[xyzinpn, ddotU, dotY] = generate_data_PLDS(T, trueparams);

% subplot(211); plot(xyzinpn.y', 'o');
% subplot(212); plot(xyzinpn.x');

%% EM algorithm

datastruct = runEMforPLDS(xyzinpn, k, pinp, p, ddotU, trueparams);

%% sanity check

nsamps = 5000;
Zn = zeros(p, T, nsamps);
Zn_estParam = zeros(p, T, nsamps);

datastruct.Mstep.inpn = trueparams.inpn;

for i=1:nsamps
    [i nsamps]
    
    sanitycheck = generate_data_PLDS(T,datastruct.Mstep);
    Zn_estParam(:,:,i) = sanitycheck.z;
    
    sanitychecktrue = generate_data_PLDS(T,trueparams);
    Zn(:,:,i) = sanitychecktrue.z;
end

%%

meanZn = mean(mean(Zn,3),2);
meanZnest = mean(mean(Zn_estParam,3),2);
 
[meanZn meanZnest]

Znrshp = reshape(Zn, p, []);
Znestrshp = reshape(Zn_estParam, p, []);

% [mean(meanZn,2)
[cov(Znrshp') cov(Znestrshp')]

%%
figure;
% subplot(2,2,[1 2]); plot([mean(mean(Zn,3),2) mean(mean(Zn_estParam,3),2)]); legend('mean of Y using true model parameters', 'mean of Y using model estimates'); xlabel('neurons'); 

% Yn_estParam = generate_data_PLDS(T,Anew,B,Cnew,d,Q,x0,V0,inpn);

% subplot(2,4,[1 2 3 4]); plot([meanZn meanZnest]); legend('mean of Y using true params', 'mean of Y using estimated params'); xlabel('neurons'); 
% set(gca,'ylim', [-3 -1.9]);

% subplot(3,3,[4 5]); plot(1:T, trueparams.x', 'k', 1:T, datastruct.Estep.mumarg', 'r');xlabel('mu of x')

% subplot(2,4,1); hinton(trueparams.A,['true A  ' num2str(max(max(trueparams.A)),'%.4f')],'standard');
% subplot(2,4,2); hinton(datastruct.Mstep.A,['Aest  ' num2str(max(max(datastruct.Mstep.A)),'%.4f')],'standard');
% subplot(2,4,3); hinton(trueparams.C,['true C  ' num2str(max(max(trueparams.C)),'%.4f')],'standard');
% subplot(2,4,4); hinton(datastruct.Mstep.C,['Cest  ' num2str(max(max(datastruct.Mstep.C)),'%.4f')],'standard');

% subplot(2,4,5); hinton(cov(Znrshp'),['cov z' num2str(max(max(cov(Znrshp'))),'%.4f')],'standard');
% subplot(2,4,6); hinton(cov(Znestrshp'),['cov z est  ' num2str(max(max(cov(Znestrshp'))),'%.4f')],'standard');

Znrshp2 = Znrshp(:,[2:end,1]);
Znestrshp2 = Znestrshp(:, [2:end,1]);

% subplot(2,4,7); hinton(cov(Znrshp', Znrshp2'),['cov(z_t, z_{t-1})' num2str(max(max(cov(Znrshp', Znrshp2'))),'%.4f')],'standard');
% subplot(2,4,8); hinton( cov(Znestrshp', Znestrshp2'),['cov(zp_t, zp_{t-1})' num2str(max(max(cov(Znestrshp', Znestrshp2'))),'%.4f')],'standard');


%% this was for comparing mine with Lars's code...

% load m;
% load c;
% load Cest;

% subplot(3,3,1); hinton(C,['true C  ' num2str(max(max(C)),'%.4f')],'standard');
% subplot(3,3,2); hinton(Cest.cmj,['Cest(mijung)  ' num2str(max(max(Cest.cmj)),'%.4f')],'standard');
% subplot(3,3,3); hinton(Cest.clr,['Cest(lars) ' num2str(max(max(Cest.clr)),'%.4f')],'standard');

% subplot(3,3,[4 5 6]); plot(1:p, m.meantrue, 'k', 1:p, m.meanmj, 'g', 1:p, m.meanlr, 'r');
% set(gca, 'ylim', [-2.1 -1.9]); legend('mean of true z', 'mijung', 'lars');

% subplot(3,3,7); hinton(c.covtrue,['cov(z true)  ' num2str(max(max(c.covtrue)),'%.4f')],'standard');
% subplot(3,3,8); hinton(c.covmj,['cov(z mijung) ' num2str(max(max(c.covmj)),'%.4f')],'standard');
% subplot(3,3,9); hinton(c.covlr,['cov(z lars)  ' num2str(max(max(c.covlr)),'%.4f')],'standard');

% subplot(3,3,7); hinton(cov(Znrshp'),['cov y' num2str(max(max(cov(Znrshp'))),'%.4f')],'standard');
% subplot(3,3,8); hinton(cov(Znestrshp'),['cov y est  ' num2str(max(max(cov(Znestrshp'))),'%.4f')],'standard');
