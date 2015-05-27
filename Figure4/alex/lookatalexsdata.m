% first session

clear all;
% close all;
clc;

session = 2; 
data = preprocessor_alex(session);

orientataion = 0;
binsize = 50;
[stims, resps] = data.get_data_for_ori(orientataion,binsize);

[p,T,r] = size(resps);


mfr = zeros(p,r);
vfr = zeros(p,r);
yy = [];

for i=1:r 
   mfr(:,i) = mean(resps(:,:,i),2);
%    vfr(:,i) = var(resps(:,:,i), [],2);
   yy = [yy resps(:,:,i)];
end

figure(1);
% subplot(211); plot(mfr');
% subplot(212); plot(var(mfr,[],2)', '.');

load signi_fallers_sort_by_varexp;
load signi_risers_sort_by_varexp;

% signi_risers_sort_by_varexp
subplot(211); plot(mfr(signi_risers_sort_by_varexp,:)')
subplot(212); plot(mfr(signi_fallers_sort_by_varexp,:)')

% figure(2);
% subplot(211); plot(mfr, '.');
% subplot(212); plot(var(yy,[],2)', '.');
% 
%%
% extract parts that have inputs
% range = floor(T/4)+1:floor(3*T/4); 
% % T = length(range);
% 
% stims = stims(:,range,1);
% resps = resps(:,range,:);
% 
% [p, T, r] = size(resps);
% 
% % total covariance
% yy = [];
% for i=1:r
%     yy = [yy resps(:,:,i)];
% end
% 
% mfr = zeros(p,r);
% vfr = zeros(p,r);
% 
% for i=1:r 
%    mfr(:,i) = mean(resps(:,:,i),2);
%    vfr(:,i) = var(resps(:,:,i), [],2);
% end


% figure(1);
% subplot(211); plot(mfr', '.');
% % subplot(212); plot(var(mfr,[],2)', '.');
% 
% 
% figure(2);
% subplot(211); plot(mfr, '.');
% subplot(212); plot(var(yy,[],2)', '.');

% subplot(211); plot(mfr(signi_risers_sort_by_varexp,:)')
% subplot(212); plot(mfr(signi_fallers_sort_by_varexp,:)')

%%

% conditional covariance
condicov = zeros(p,p,r);
for i=1:r
    condicov(:,:,i) = cov(resps(:,:,i)',1);
end
% subplot(2,2,[1 2]); imagesc([cov(yy') mean(condicov,3)]); colorbar

meanmat = zeros(p,r);
for i=1:r
    meanmat(:,i) = mean(resps(:,:,i),2);
end
condicov2 = cov(meanmat',1);


[ sum(sum(cov(yy'))) sum(sum(mean(condicov,3))) sum(sum(condicov2))] 

%    10.7755    8.5244    2.2490
sum(sum((cov(yy') - (mean(condicov,3) + condicov2)).^2))

%%

numtimebins = T;
autocorr = zeros(p, numtimebins+1);

for whichcell = 1:p
    autocorr(whichcell, :) = xcov(yy(whichcell,:), numtimebins/2, 'unbiased');
end

avgautocorr_acrcells = mean(autocorr);

% figure; 
subplot(211); plot(autocorr'); title('total auto-cov of each cell');
subplot(212); plot(avgautocorr_acrcells); title('sum of total auto-cov');

%%

autocorr_condi = zeros(p, numtimebins+1, r);

for whichcell = 1:p
    
    for whichtrial = 1:r
        autocorr_condi(whichcell, :, whichtrial) = xcov(resps(whichcell,:,whichtrial), numtimebins/2, 'unbiased');
    end
    
end

%%
autocorr_per_eachtrial = squeeze(mean(autocorr_condi));
avgautocorr_condi = mean(autocorr_per_eachtrial,2);

% figure;
plot(211);
plot(1:numtimebins+1, avgautocorr_acrcells/max(avgautocorr_acrcells), 'k', 1:numtimebins+1, avgautocorr_condi/max(avgautocorr_condi), 'r--')
legend('total covariance', 'conditional covariance');
%% 

subplot(211);
plot(mean(yy'), 'o')

%%

% save alexdata_session3_org135_bin40.mat


