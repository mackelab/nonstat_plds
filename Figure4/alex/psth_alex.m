function psth_alex(bin_size)

n_data_sets = preprocessor_alex.available_sessions();
for data_set = 1:n_data_sets
    data = preprocessor_alex(data_set);
    
    % Get data and dimensions
    x = bin_size:bin_size:(2 * data.stim_duration);
    y = 1000 * data.get_binned_spikes(bin_size, false) / bin_size;
    
    [C, T, N] = size(y);
    
    num_cond = length(data.directions);
    assert(num_cond == 16, 'Not implemented, yet!');
    
    % Compute PSTH for each cell and each orientation
    psth = zeros(C, T, num_cond);
    for cond = 1:num_cond
        psth(:, :, cond) = sum(y(:, :, data.conditions == cond), 3) / sum(data.conditions == cond);
    end
    
    % Plot all
    figure('Name', sprintf('Dataset %d - All', data_set))
    p = panel();
    p.pack(4, 4);
    
    max_y = max(max(max(psth)));
    for cond = 1:num_cond
        
        [a, b] = ind2sub([4, 4], cond);
        p(a, b).select()
        
        plot(x, psth(:, :, cond)');
        ylim([0, max_y]);
        ylabel('spikes / sec')
        xlim([0, x(end) + 1]);
        xlabel('time (ms)');
        title(sprintf('%.1f degree', 180 * data.directions(cond) / pi));
    end
    
    mu = squeeze(mean(psth, 1));
    sem = squeeze(std(psth, 1)) / sqrt(size(psth, 1));
    
    % Plot mean
    figure('Name', sprintf('Dataset %d - Mean', data_set))
    p = panel();
    p.pack(4, 4);
    
    min_y = min(min(mu - sem));
    max_y = max(max(mu + sem));
    for cond = 1:num_cond
        
        [a, b] = ind2sub([4, 4], cond);
        p(a, b).select()
        
        errorbar(x, mu(:, cond), sem(:, cond));
        ylim([min_y, max_y]);
        ylabel('spikes / sec')
        xlim([0, x(end) + 1]);
        xlabel('time (ms)');
        title(sprintf('%.1f degree', 180 * data.directions(cond) / pi));
    end
    
    
end