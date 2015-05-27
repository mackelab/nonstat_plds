% first session

% clear all;
% close all;
% clc;
addpath ../../../gpml-matlab/gpml/
addpath ../alex/

%work at binsize 50ms, and ignore pre and post stimulus data:

binsize = 50;

%load all three sessions
for session=1:4
    session
    %collect session
    data(session) = preprocessor_alex(session);
    
    %loop over all directions
    for dir_index=1:16
        dir_index
        direction=round(data(session).directions(dir_index)/2/pi*360*2)/2;
        [stims, resps] = data(session).get_data_for_ori(direction,binsize);
        n_neurons=size(resps,1);
        
        rates=squeeze(mean(resps,2));
        for i=1:n_neurons
            %for each neuron and each orientation, we do four things
            
            stats(session).mean(i,dir_index)=mean(rates(i,:));
            stats(session).var(i,dir_index)=var(rates(i,:),0);
            stats(session).rates(i,dir_index,:)=rates(i,:);
            %a) we do a t-test on whether the firing rate on the first 50 trials is
            %different to the second 50 trials
            
            [junk, stats(session).mean50.p_value(i,dir_index)]=ttest2(rates(i,1:50),rates(i,51:end),.05,'both');
            stats(session).mean50.means(i,dir_index,1)=mean(rates(i,1:50));
            stats(session).mean50.means(i,dir_index,2)=mean(rates(i,51:end));
            stats(session).mean50.vars(i,dir_index,1)=var(rates(i,1:50));
            stats(session).mean50.vars(i,dir_index,2)=var(rates(i,51:end));
            stats(session).mean50.dprime(i,dir_index)=diff(stats(session).mean50.means(i,dir_index,:))/sqrt(mean(stats(session).mean50.vars(i,dir_index,:)));
           
        
            %b) we fit a linear regression
            x=[ones(100,1),[1:100]'/100];
            [B,BINT,R,RINT,STATS] = regress(rates(i,:)',x);
            stats(session).linreg.coeffs(i,dir_index,:)=B;
            stats(session).linreg.p_value(i,dir_index)=STATS(3);
            stats(session).linreg.var_exp(i,dir_index)=1-var(R,0)/var(rates(i,:),0);
            if stats(session).linreg.var_exp(i,dir_index)<-.01;
                keyboard
            end
            
            %c) we fit a quadratic regression (without a linear part)
            x=[ones(100,1),([1:100]'/100).^2];
            [B,BINT,R,RINT,STATS] = REGRESS(rates(i,:)',x);
            stats(session).quadreg.coeffs(i,dir_index,:)=B;
            stats(session).quadreg.p_value(i,dir_index)=STATS(3);
            stats(session).quadreg.var_exp(i,dir_index)=1-var(R,0)/var(rates(i,:),0);
            %In each case, we store both p-values and (some crude) measure of effect
            %sizes.
            
            %d) full quadratic regression
            x=[ones(100,1),([1:100]'/100),([1:100]'/100).^2];
            [B,BINT,R,RINT,STATS] = REGRESS(rates(i,:)',x);
            stats(session).fullquadreg.coeffs(i,dir_index,:)=B;
            stats(session).fullquadreg.p_value(i,dir_index)=STATS(3);
            stats(session).fullquadreg.var_exp(i,dir_index)=1-var(R,0)/var(rates(i,:),0);
            
            %now, lets use GP regression to find the best smoothing kernel for each
            %neuron. This is probably a bit of an overkill but maybe important and
            %useful later on:
            x=[[1:100]'];
            y=rates(i,:)'-mean(rates(i,:)');
            covfunc = {'covSum', {'covSEiso','covNoise'}};
            loghyper = [log(1.0); log(1.0); log(0.1)];
            loghyper = minimize(loghyper, 'gpr', -100, covfunc, x, y);
            stats(session).gp.tau(i,dir_index)=exp(loghyper(1));
            stats(session).gp.var(i,dir_index)=exp(loghyper(2));
            stats(session).gp.noisevar(i,dir_index)=exp(loghyper(3));
            [mu] = gpr(loghyper, covfunc, x, y,x);
            stats(session).gp.smoothrate(i,dir_index,:)=mu+mean(rates(i,:)');
            %keyboard
            
        end
    end
end

save stats stats
%%
clear all, close all
load stats
h(1)=figure;

cutoff=0.05
%display significant p-values for each session

for session=1:4
    for i=1:3
    subplot(3,4,session+4*(i-1))
    switch i
        case 1
            d=stats(session).mean50.p_value<cutoff;
        case 2
            d=stats(session).linreg.p_value<cutoff;
        case 3
            d=stats(session).quadreg.p_value<cutoff;
    end
    imagesc(d);
    title(num2str(mean(d(:))));
    end
end


%% from these data, it very clearly looks like we want to look at session 2,
%so lets look at this in more detail. 
%It also looks as if the effects are similar across orientations, so lets
%just take the first one:
session=2;
dir_index=1;


h(2)=figure;


%sort neurons by variance explained: 
[varexp_sorted,sort_by_varexp]=sort(stats(session).linreg.var_exp(:,dir_index),'descend');

subplot(3,3,1)
select=stats(session).linreg.p_value(:,dir_index)<cutoff & stats(session).linreg.coeffs(:,dir_index,2)>0;
d=squeeze(stats(session).rates(:,dir_index,:));
d(~select,:)=0;
imagesc(d(sort_by_varexp,:))
signi_risers=find(select);
signi_risers_sort_by_varexp=sort_by_varexp(select(sort_by_varexp));


title('FRs of signi risers')

subplot(3,3,2)
select=stats(session).linreg.p_value(:,dir_index)<cutoff & stats(session).linreg.coeffs(:,dir_index,2)<0;
d=squeeze(stats(session).rates(:,dir_index,:));
d(~select,:)=0;
imagesc(d(sort_by_varexp,:))
signi_fallers=find(select);
signi_fallers_sort_by_varexp=sort_by_varexp(select(sort_by_varexp));

title('FRs of signi fallers')

subplot(3,3,3)
select=stats(session).linreg.p_value(:,dir_index)>cutoff ;
d=squeeze(stats(session).rates(:,dir_index,:));
d(~select,:)=0;
imagesc(d(sort_by_varexp,:))
signi_non=find(select);
signi_non_sort_by_varexp=sort_by_varexp(select(sort_by_varexp));
title('FRs of non-sign neurons')


subplot(3,3,4)
d=squeeze(stats(session).rates(signi_risers_sort_by_varexp,dir_index,:));
plot(d(1:min(end,10),:)');
title('Firing rate of top 10 risers')

subplot(3,3,5)
d=squeeze(stats(session).rates(signi_fallers_sort_by_varexp,dir_index,:));
plot(d(1:min(end,10),:)');
title('Firing rate of top 10 fallers')

subplot(3,3,6)
hist(varexp_sorted,20)
title('Histogram of variance explained')

subplot(3,3,7)
d=squeeze(stats(session).gp.smoothrate(signi_risers_sort_by_varexp,dir_index,:));
plot(d(1:min(end,10),:)');
title('Smoothed Firing rate of top 5 risers')

subplot(3,3,8)
d=squeeze(stats(session).gp.smoothrate(signi_fallers_sort_by_varexp,dir_index,:));
plot(d(1:min(end,10),:)');
title('Smoothed rate of top 5 fallers')


save stats stats signi* varexp*






