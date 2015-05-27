function barplot_alex()
    n_data_sets = preprocessor_alex.available_sessions();
    hold on;
    for data_set = 1:n_data_sets
        data = preprocessor_alex(data_set);
        trial_length = (2 * data.stim_duration + 1000) / 1000; % in s
        spikes_per_s = sum(cellfun(@length, data.spikes(:, :)), 1) / (size(data.spikes, 1) * trial_length);
        boxplot(spikes_per_s, 'positions', data_set, 'widths', 0.8);
    end
    hold off
    xlim([0, n_data_sets + 1])
    set(gca, 'XTick', 1:n_data_sets)
    set(gca, 'XTickLabel', 1:n_data_sets)
    xlabel('session index')
    ylabel('spikes / sec')