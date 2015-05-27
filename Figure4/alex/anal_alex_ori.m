function [data, plot_handle] = anal_alex_ori(varargin)

seed = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(seed);

%% Config
% ToDo: Make configurable
N_RANDITER = 24;
N_ITER = 100;


inp_parse = inputParser();

inp_parse.addRequired('recording_index');
inp_parse.addRequired('orientation');
inp_parse.addRequired('bin_size');

inp_parse.addOptional('run_ml', false);
inp_parse.addOptional('width_hist', 20);
inp_parse.addOptional('n_filter_features', 5);
inp_parse.addOptional('a0', 1e-2);
inp_parse.addOptional('b0', 1e-4);
inp_parse.addOptional('N_RANDSTART', N_RANDITER);
inp_parse.addOptional('N_ITER', N_ITER);

inp_parse.parse(varargin{:});

struct2ws('inp_parse.Results');

%% Load data
pp = preprocessor_alex(recording_index);
[stim_x, y] = pp.get_data_for_ori(orientation, bin_size);

[D_stim, T, N] = size(stim_x);
C = size(y, 1);


%% Generate constant features
constant_x = ones(1, T, N);

full_x = [constant_x; stim_x];

state_x = arrayfun(@(~) full_x, 1:C, 'UniformOutput', false);

%% Generate filters (by referencing to the same filter)
hist_filter = arrayfun(@(~) filter_postspike(width_hist, n_filter_features), 1:C, 'UniformOutput', false);

%% Create states
default_state_prior = prior_gauss(prior_gamma(a0 * ones(n_filter_features + 1, 1), ...
    b0 * ones(n_filter_features + 1, 1)));
default_state_predictor = predictor_logreg(default_state_prior);

pred_up_high = arrayfun(@(~) default_state_predictor.copy(), 1:C, 'UniformOutput', false);
s_up_high = node_state('s_up_high', [pred_up_high{:}]);

pred_up_low = arrayfun(@(~) default_state_predictor.copy(), 1:C, 'UniformOutput', false);
s_up_low = node_state('s_up_low', [pred_up_low{:}]);

pred_down = arrayfun(@(~) default_state_predictor.copy(), 1:C, 'UniformOutput', false);
s_down = node_state('s_down', [pred_down{:}]);

%% Create nodes
default_gate_prior = prior_gauss(prior_gamma(a0 * ones(1 + 3, 1), b0 * ones(1 + 3, 1)));
default_gate_predictor = predictor_logreg(default_gate_prior);

gate_high_low = node_gate_hmm('g_high_low', default_gate_predictor.copy(), s_up_high, s_up_low);

gate_up_down = node_gate_hmm('g_up_down', default_gate_predictor.copy(), gate_high_low, s_down);

%% Initialize data struct

data = node_data.init(gate_up_down, y, ...
    'g_up_down', constant_x, [], ...
    'g_high_low', constant_x, [], ...
    's_up_high', state_x, hist_filter,...
    's_up_low', state_x, hist_filter, ...
    's_down', state_x, hist_filter);

%% Fit VB
resp_state_start = cell(1, N_RANDSTART);
resp_joined_start = cell(1, N_RANDSTART);

bayes_loglike = nan(N_RANDSTART,N_ITER);
bayes_loglike_point = nan(N_RANDSTART,N_ITER);
bayes_evidence = nan(N_RANDSTART,N_ITER);
bayes_ud = arrayfun(@(~) gate_up_down.copy, 1:N_RANDSTART, 'UniformOutput', false);

for n = 1:N_RANDSTART
    display(['Bayes: fitting random itialization ' num2str(n)])
    resp_joined_start{n} = rand(3,N*(T-1),3) / (3^2);
    resp_state_start{n} = rand(3,N*T);
    temp_model = bayes_ud{n};
    [resp_state, resp_joined, logexpect] = ...
        temp_model.update(data, resp_state_start{n}, resp_joined_start{n});
    for i = 1:N_ITER-1
        display(['Iteration ', num2str(i)])
        [resp_state, resp_joined, logexpect] = ...
            temp_model.update(data, resp_state, resp_joined);
        bayes_loglike(n,i) = sum(logexpect);
        bayes_evidence(n, i) = temp_model.evidence_bound(data, bayes_loglike(n,i));
        bayes_loglike_point(n, i) = sum(temp_model.loglike(data));
    end
    bayes_ud{n} = temp_model;
end
[~, best_bayes] = max(bayes_evidence(:,end-1));
bayes_up_down = bayes_ud{best_bayes};

bayes_path = reshape(bayes_up_down.decode_state(data), 3, T, N);

%% Fit ML

% if(run_ml)
%     ml_up_down = gate_up_down.copy();
%     arrayfun(@(p) set(p, 'prior', prior_ml(zeros(size(p.prior.mean)))), ...
%         [ml_up_down.get_nodes.predictor])
%     
%     
%     ml_loglike = nan(N_RANDSTART,N_ITER);
%     ml_loglike_point = nan(N_RANDSTART,N_ITER);
%     ml_evidence = nan(N_RANDSTART,N_ITER);
%     ml_ud = arrayfun(@(~) ml_up_down.copy, 1:N_RANDSTART, 'UniformOutput', false);
%     
%     parfor n = 1:N_RANDSTART
%         display(['ML: fitting random itialization ' num2str(n)])
%         temp_model = ml_ud{n};
%         [resp_state, resp_joined, logexpect] = ...
%             temp_model.update(data, resp_state_start{n}, resp_joined_start{n});
%         for i = 1:N_ITER-1
%             [resp_state, resp_joined, logexpect] = ...
%                 temp_model.update(data, resp_state, resp_joined);
%             ml_loglike(n,i) = sum(logexpect);
%             ml_evidence(n, i) = temp_model.evidence_bound(data, ml_loglike(n,i));
%             ml_loglike_point(n, i) = sum(temp_model.loglike(data));
%         end
%         ml_ud{n} = temp_model;
%     end
%     [~, best_ml] = max(ml_evidence(:,end-1));
%     ml_up_down = ml_ud{best_ml};
%     
%     ml_path = reshape(ml_up_down.decode_state(data), 3, T, N);
% end

%% Save results
save(sprintf('alex%d_ori%.1f_bs%d.mat', recording_index, orientation, bin_size));

data = ws2struct();
plot_handle = @plot_data;

end


function plot_data(data)

struct2ws('data');

%% Plot data
% color_low   = [0.8000    0.9216    0.7725];
% color_medium = [251 154 153] ./ 255;
% color_high = [255 26 28]./255;

color_low   = [77,175,74] ./ 255;
color_medium = [255,127,0] ./ 255;
color_high = [228,26,28]./255;

color_all = [color_high; color_medium; color_low];

[fig_pars, axe_pars] = plot_pub_pars();

figure(fig_pars);
trial_to_plot = 1; %randi(N);

p = panel();
p.pack('v', {1/7 []})
p(2).pack(5)

% Plot stimulus
p(1).select();

hold on
plot(squeeze(stim_x(:, :, trial_to_plot)'), 'k');
ylim([-1.2, 1.2])
set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'YTick', [])
ylabel('Stimulus')
box off


for trial_to_plot = 1:5
% Plot states and spikes (bayes)
p(2,trial_to_plot).select();

plot_states(logical(bayes_path(:, :, trial_to_plot)), {color_high, color_medium, color_low}, ...
    {'up - high', 'up - low', 'down'});

hold on;
height = 1 / C;
for t = 1:T
    for i = 1:C
        if y(i, t, trial_to_plot)
            plot([t, t], [(i - 1) * height, i * height], 'k', 'LineWidth', 1)
        end
    end
end
set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'YTick', [], axe_pars)
ylabel('Bayes')

plot(1:T, mean(y(:,:,trial_to_plot),1), 'color', [0.4 0.4 0.4], 'LineWidth', 2)

% Plot states and spikes (ML)
if(run_ml)
    p(3).select();
    
    plot_states(logical(ml_path(:, :, trial_to_plot)), {color_high, color_medium, color_low}, ...
        {'up - high', 'up - low', 'down'});
    
    hold on;
    height = 1 / C;
    for t = 1:T
        for i = 1:C
            if y(i, t, trial_to_plot)
                plot([t, t], [(i - 1) * height, i * height], 'k', 'LineWidth', 1)
            end
        end
    end
    set(gca, 'XTickLabel', {}, 'YTickLabel', {}, 'YTick', [], axe_parsaxe_pars)
    ylabel('ML')
end
end

%%
figure;
plot([mean(mean(y(:, bayes_path(1,:,:) == 1))), ...
    mean(mean(y(:, bayes_path(2,:,:) == 1))), ...
    mean(mean(y(:, bayes_path(3,:,:) == 1)))])

%%
data_autocov = autocovariance(y(:, T/4+1:3/4*T, :), bayes_path(:, T/4+1:3/4*T, :), 50);
data_autocov_innocent = autocovariance(y(:, T/4+1:3/4*T, :), ones(1,T/2,N), 50);

ylim_use = 1.1*[min(min(data_autocov(2:end,:))), max(max(data_autocov(2:end,:)))];
figure
for c = 1:C
    subplot(ceil(sqrt(C)),ceil(sqrt(C)), c)
    
    hold all
    for i = 1:3
        plot(1:size(data_autocov,1)-1, data_autocov(2:end,i,c), ...
            'Color', color_all(i,:), 'LineWidth', 2)
    end
    plot(1:size(data_autocov,1)-1, data_autocov_innocent(2:end,1,c), ...
            'k--', 'LineWidth', 2)
    axis off
    ylim(ylim_use)
end

end