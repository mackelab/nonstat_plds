clear all;
close all;
clc;

addpath ../../gpml-matlab/gpml/

%work at binsize 50ms, and ignore pre and post stimulus data:

binsize = 50;

session=2;

% load data for the session
data = preprocessor_alex(session);

% think about one direction
dir_index=1;

direction=round(data.directions(dir_index)/2/pi*360*2)/2;
[stims, resps] = data.get_data_for_ori(direction,binsize);
n_neurons=size(resps,1);

%%
rates=squeeze(mean(resps,2));
for i=1:n_neurons
    %for each neuron and each orientation, we do four things
    
    stats.mean(i,dir_index)=mean(rates(i,:));
    stats.var(i,dir_index)=var(rates(i,:),0);
    stats.rates(i,dir_index,:)=rates(i,:);
    
    %now, lets use GP regression to find the best smoothing kernel for each
    %neuron. This is probably a bit of an overkill but maybe important and
    %useful later on:
    x=[[1:100]'];
    y=rates(i,:)' - mean(rates(i,:)');
    covfunc = {'covSum', {'covSEiso','covNoise'}};
    loghyper = [log(1.0); log(1.0); log(0.1)];
    loghyper = minimize(loghyper, 'gpr', -100, covfunc, x, y);
    stats.gp.tau(i,dir_index)=exp(loghyper(1)); % this is exp(-(xi-xj)^2/tau), length scale
    stats.gp.var(i,dir_index)=exp(loghyper(2));
    stats.gp.noisevar(i,dir_index)=exp(loghyper(3));
    [mu] = gpr(loghyper, covfunc, x, y,x);
    stats.gp.smoothrate(i,dir_index,:)= mu + mean(rates(i,:)');
    %keyboard
    
end

%%


