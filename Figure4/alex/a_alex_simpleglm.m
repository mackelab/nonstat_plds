clearvars

seed = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(seed);

%% Config
run_ml = false;
logexp = @(x) 1./(1 + exp(-x));
%% Load data
bin_size = 10; % in ms
recording_index = 2;
orientation = 90;

pp = preprocessor_alex(recording_index);
[stim_x, y] = pp.get_data_for_ori(orientation, bin_size);

[D_stim, T, N] = size(stim_x);
C = size(y, 1);


%% Generate constant features
state_x = arrayfun(@(~) stim_x, 1:C, 'UniformOutput', false);

%% Generate filters (by referencing to the same filter)
width_hist = 20;
n_filter_features = 5;
hist_filter = arrayfun(@(~) filter_postspike(width_hist, n_filter_features), 1:C, 'UniformOutput', false);

%% Create states
a0 = 1e-2;
b0 = 1e-4;

default_state_prior = prior_gauss(prior_gamma(a0 * ones(n_filter_features + D_stim, 1), ...
    b0 * ones(n_filter_features + D_stim, 1)));
default_state_predictor = predictor_logreg(default_state_prior);

state_predictor = arrayfun(@(~) default_state_predictor.copy(), 1:C, 'UniformOutput', false);
bayes_state = node_state('GLM', [state_predictor{:}]);

%% Initialize data struct

data = node_data.init({'GLM'}, y, ...
    'GLM', state_x, hist_filter);

%% Fit VB

z = ones(1, size(y,2) * size(y,3));
for i = 1:10
    display(['Iteration ', num2str(i)])
    bayes_state.update(data,z);
    logexpect = bayes_state.expectation(data);
    bayes_loglike(i) = sum(logexpect);
    bayes_evidence(i) = bayes_state.evidence_bound() + bayes_loglike(i);
%     bayes_loglike_point(i) = sum(bayes_state.loglike(data));
end

% %%
% figure;
% plot(bayes_evidence)
% 
% figure;
% plot(bayes_loglike)
% 
% 
% figure;
% hold all
% ha = arrayfun(@(i) plot(1:width_hist, bayes_state.predictor(i).prior.mean(D_stim+1:end)' * hist_filter{i}.filter, ...
%     'LineWidth', 2), 1:length(bayes_state.predictor));
% 
% 
% figure
% hold all
% ha = arrayfun(@(i) plot(logexp(bayes_state.predictor(i).prior.mean(1:D_stim)' * stim_x(:,:,1)), ...
%     'LineWidth', 2), 1:length(bayes_state.predictor));
% 
% 
% %%
% p_hist = mean(y,3);
% 
% ylim_use = [0, 1.1 * max(max(p_hist))];
% 
% gfilter = fspecial('gaussian',[1 5], 1);
% p_diff = zeros(C, 2);
% figure;
% stim_id = find(stim_x(1,:,1));
% for c = 1:C
% %     [i, j] = ind2sub([ceil(sqrt(C)), ceil(sqrt(C))], c);
% %     p(j, i).select();
%     subplot(ceil(sqrt(C)),ceil(sqrt(C)), c)
%     
%     hold all
%     p_bin = p_hist(c,:);
%     
%     p_pred = mean(reshape(bayes_state.predictor(c).predict(data.trial_x{c}), T, []), 2)';
%     
% %     p_pred = logexp(bayes_state.predictor(c).prior.mean(1:D_stim)' * stim_x(:,:,1));
%     p_smooth = imfilter(p_hist(c,:), gfilter, 'replicate');
%     bH = bar(p_bin);
%     h = findobj(bH,'Type','patch');
%     set(h,'FaceColor','k','EdgeColor','k','facealpha',0.25, 'edgealpha', 0)
%     plot(p_smooth, 'LineWidth', 2)
%     plot(p_pred, 'LineWidth', 2)
%     p_diff(c, :) = [(mean(p_pred(stim_id)) - mean(p_bin(stim_id))) /mean(p_bin(stim_id)), ...
%         (mean(p_smooth(stim_id)) - mean(p_bin(stim_id))) / mean(p_bin(stim_id))];
%     axis tight
%     ylim(ylim_use)
%     axis off
% end
% 
% figure; hold all
% boxplot(p_diff)
% plot(xlim, [0,0], 'k--')
% 
% figure;
% for c = 1:C
% %     [i, j] = ind2sub([ceil(sqrt(C)), ceil(sqrt(C))], c);
% %     p(j, i).select();
% subplot(ceil(sqrt(C)),ceil(sqrt(C)), c)
%     hold all
%     plot(1:width_hist, bayes_state.predictor(c).prior.mean(D_stim+1:end)' * hist_filter{c}.filter, ...
%     'LineWidth', 2)
% end