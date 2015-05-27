function corr_alex(bin_size)

n_data_sets = preprocessor_alex.available_sessions();
for data_set = 1:n_data_sets
    data = preprocessor_alex(data_set);
    
    % Get data and dimensions
    y = data.get_binned_spikes(bin_size, false);
    
    C = size(y, 1);
    
    num_cond = length(data.directions);
    assert(num_cond == 16, 'Not implemented, yet!');
    
    % Compute PSTH for each cell and each orientation
    corr = zeros(C, C, num_cond);
    for cond = 1:num_cond
        corr(:, :, cond) = corrcoef(reshape(y(:, :, data.conditions == cond), C, [])');
    end
    
    % Plot all
    figure('Name', sprintf('Dataset %d', data_set))
    colormap('hot')
    p = panel();
    p.pack(4, 4);
    
    for cond = 1:num_cond
        
        [a, b] = ind2sub([4, 4], cond);
        curr_p = p(a, b);
        curr_p.select()
        curr_p.margin = [20, 20, 0, 0];
        
        imagesc(corr(:, :, cond));
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        colorbar();
        ylim([1, C])
        xlim([1, C])
        title(sprintf('%.1f degree', 180 * data.directions(cond) / pi));
    end  
    drawnow()
end