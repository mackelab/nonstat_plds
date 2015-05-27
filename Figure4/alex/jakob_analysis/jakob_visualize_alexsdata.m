
%%
clear all, close all
load stats
h(1)=figure;
addpath ~/ME/Proj/mackelab_stuff/code/matlab/plot/
cutoff=0.05
%display significant p-values for each session

for session=1:4
    for i=1:3
    subplot(3,4,session+4*(i-1))
    switch i
        case 1
            d=stats(session).mean50.p_value<cutoff;
            imagesc(d);
            title(['Mean signi: ', num2str(mean(d(:)))]);
        case 2
            d=stats(session).linreg.p_value<cutoff;
            imagesc(d);
            title(['Lin Reg signi: ', num2str(mean(d(:)))]);
        case 3
            d=stats(session).quadreg.p_value<cutoff;
            imagesc(d);
            title(['Quad Reg signi: ', num2str(mean(d(:)))]);
    end
    
    end
end

PrintFigure(h(1),'pdf','Overview_Datasets',[],{'Papersize',[30,20],'PaperPosition',[0,0,30,20]});


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

subplot(3,3,9)

d= stats(session).gp.tau(:,dir_index);
hist(min(d,50))
title('Histogram of taus')

%one other way to get at 'nonstationarity' neurons are ones for which the
%gp-estimate of tau is not very very huge:

signi_tau_small=find(d<50);


save stats stats signi* varexp*

PrintFigure(h(2),'pdf','Dataset_2',[],{'Papersize',[30,20],'PaperPosition',[0,0,30,20]});





