%% generate data for Figure 2
% April 22, 2014
% wrote by Mijung Park

clear all;
close all;
clc;

addpath ../testcode_nonstationaryPLDS_A_multipleRecordings/

% essential quantities (fixing pinp==0)
k = 4; % dimensionality of hidden states
p = 40; % output dimension (e.g., # of neurons)
pinp = 2; % input dimension
T = 200; % time 1:T (for training)
r = 100; % # of recordings

%% generate params (A, B, C, h, d)

% (1) A
A = zeros(k); A(3,3) = 0.9; A(4,4) = 0.9;

% (2) B
B = zeros(k, pinp); B(1:2,1:2) = [5 0; 0 6]; 

% (3) C
C = zeros(p,k); 
maxC1 = 0.1;
maxC2 = 0.1;
C(1:p/2,1) = maxC1*ones(p/2,1); C(1:p/2,3) = maxC1*ones(p/2,1);
C(p/2+1:p,2) = maxC2*ones(p/2,1); C(p/2+1:p,4) = maxC2*ones(p/2,1);
 
C(:,3) = C(:,3)/3;
C(:,4) = C(:,4)/3;

% C(1:p/2,1) = maxC*randn(p/2,1); C(1:p/2,3) = maxC*randn(p/2,1);
% C(p/2+1:p,2) = maxC*randn(p/2,1); C(p/2+1:p,4) = maxC*randn(p/2,1);


% (4) d
% d = -2*ones(p,1);
d = [-2.2*ones(p/2,1); -1.2*ones(p/2,1)];

% (5) h
h = zeros(k,r);
h13 = normpdf(1:r, r/4, r/8);
h13 = 12*h13/max(h13);
h24 = - fliplr(h13);
h(1,:) = h13; h(3,:) = h13;
h(2,:) = h24; h(4,:) = h24; 

% subplot(211); plot(1:r, h13);
% subplot(212); plot(1:r, h24);

% (6) priors on x and noise variance for x
x0 = zeros(k, 1); V0 = eye(k); sqQ = sqrtm(eye(k));

% generate inputs

twoperiods = 100;
% inp - u by Tn input sequence
ti = [1:T]/twoperiods*4*pi;

% I don't know why he does it this way, but
% apparently this is how he generated inputs (u_t for each t)
window = [zeros(1,T/4) ones(1,T/2) zeros(1,T/4)];
inpn = 0.4*[cos(ti).*window; sin(ti).*window]; % assumes u = 2;

% subplot(211); plot(inpn(1,:))
% subplot(212); plot(inpn(2,:))
% % putting parameters to params structure
params.A = A;
params.B = B;
params.C = C;
params.d = d;
params.Q = sqQ;
params.h = h;
params.inpn = inpn;
params.x0 = x0;
params.V0 = V0;

Model = 'NSFR';
xyzinpn = generate_data_PLDS_multiple_recordings(T,params,r,Model);

meanspikecount = zeros(p, r);
meanz = zeros(p,r);
meanx = zeros(k,r);
y = zeros(p, T, r);
z = zeros(p, T, r);

for i=1:r
    y(:,:,i) = xyzinpn{i}.y;  
    x = xyzinpn{i}.x;
    meanspikecount(:,i) = mean(xyzinpn{i}.y,2);
    meanz(:,i) = mean(xyzinpn{i}.z,2);
    meanx(:,i) = mean(x,2);
    z(:,:,i) = xyzinpn{i}.z;
end


% save all_NSFR.mat

% load all_NSFR.mat;

% subplot(311); plot(1:r, meanz(1:p/2,:)', 'r', 1:r, meanz(p/2+1:p,:)', 'b'); ylabel('mean of z')
% subplot(312); plot(1:r, meanspikecount(1:p/2,:)', 'r', 1:r, meanspikecount(p/2+1:p,:)', 'b'); ylabel('mean of y')
% cov

% total covariance
yy = [];
zz = [];
for i=1:r
    yy = [yy y(:,:,i)];
    zz = [zz z(:,:,i)];
end

% subplot(313); plot([mean(yy(1:p/2,:)') mean(yy(1+p/2:p,:)')])

% conditional covariance
condicov = zeros(p,p,r);
for i=1:r
    condicov(:,:,i) = cov(y(:,:,i)',1);
end
% subplot(2,2,[1 2]); imagesc([cov(yy') mean(condicov,3)]); colorbar

meanmat = zeros(p,r);
for i=1:r
    meanmat(:,i) = mean(y(:,:,i),2);
end
condicov2 = cov(meanmat',1);


% subplot(223); 
% figure;
% plot(mean(condicov,3), cov(yy'),'bo'), line([0,0.5],[0,0.5]); axis image
% xlabel('conditional cov'); ylabel('total covariance');

    
sum(sum((cov(yy') - (mean(condicov,3) + condicov2)).^2))

%% saved names

% save all_NSFR.mat
% save all_NSFR_lowerRates.mat

%%
% plot(y(:,:,1)', '.')
% binsize = max(max(max(y)));
% binsize = 500; % ms

% to fix spikerate to 10Hz
% 1s : 10 spikes = binsize : mean(mean(meanspikecount))
% binsize = 0.1*mean(mean(meanspikecount))

% in miliseconds
% binsize =  round(0.2*mean(mean(meanspikecount))*1000)
binsize = 50
% to make spike rate 5Hz

rasterplot = zeros(p, binsize*T, r); 
    
for trial = 1:r   
    for whichneuron = 1: p      
        for t=1:T
            
%             [trial whichneuron t]
            
            if y(whichneuron,t,trial)>0
                howmanyspikes = y(whichneuron,t,trial);
                randomspikelocations = rand(howmanyspikes,1);
%                 middle_bin = (t-1)*binsize + binsize/2;
                spikelocations_ineachbin = 1 + (t-1)*binsize + floor(binsize.*randomspikelocations);
                
                % how many spikes?
                rasterplot(whichneuron, spikelocations_ineachbin, trial) = 1;                
            end
        end
    end
end


%%
figure(2); subplot(2,3,[ 4 5 6]);

hold on;
% Plot spikes
height = 1 / p;
trial_to_plot = 75; 

for i=1:p
    for t = 1:binsize*T
        if rasterplot(i, t, trial_to_plot)
            plot([t, t], [(i - 1) * height, i * height], 'k', 'LineWidth', 1); 
%             ylabel('corr= -0.1, trial# 50');
        end
    end
end
hold off

%%

% addpath ../mattBealsCode_v3_4_1/
% subplot(221); hinton(params.A,['true A  ' num2str(max(max(params.A)),'%.4f')],'standard');
% subplot(222); hinton(params.B,['true B  ' num2str(max(max(params.B)),'%.4f')],'standard');
% subplot(223); hinton(params.C,['true C  ' num2str(max(max(params.C)),'%.4f')],'standard');

%%

% numtimebins = size(xcorr(yy(1,:)), 2);
numtimebins = T/4;
autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
%      autocorr(whichcell, :) = xcov(yy(whichcell,:), numtimebins/2, 'unbiased');
    autocorr(whichcell, :) = xcov(zz(whichcell,:), numtimebins/2, 'biased');
end

avgautocorr_acrcells = mean(autocorr);

% figure; plot(autocorr'); title('total auto-cov of each cell');

autocorr_condi = zeros(p, numtimebins+1, r);

for whichcell = 1:p
    
    for whichtrial = 1:r
        autocorr_condi(whichcell, :, whichtrial) = xcov(z(whichcell,:,whichtrial), numtimebins/2, 'biased');
%         autocorr_condi(whichcell, :, whichtrial) = xcov(y(whichcell,:,whichtrial), numtimebins/2, 'biased');
    end
    
end

autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

figure
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r--')
legend('total covariance', 'conditional covariance');
%% 




















