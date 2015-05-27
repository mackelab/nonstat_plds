function binsize_alex(data_set)

bin_sizes = [2:5, 10:10:50, 100];

mean_lost = nan(length(bin_sizes), 1);
sem_lost = nan(length(bin_sizes), 1);

data = preprocessor_alex(data_set);

n_spikes = sum(cellfun(@length, data.spikes(:, :)), 1);
for i = 1:length(bin_sizes)
    display(['Counting spike for a bin size of ', num2str(bin_sizes(i)), '...']);
    y = data.get_binned_spikes(bin_sizes(i), false);
    
    percentage_per_cell = nan(length(n_spikes), 1);
    for unit = 1:length(n_spikes)
        percentage_per_cell(unit) = 100 * sum(y(unit, y(unit, :, :) > 1) - 1) / n_spikes(unit);
    end
    
    mean_lost(i) = mean(percentage_per_cell);
    sem_lost(i) = std(percentage_per_cell) / sqrt(length(percentage_per_cell));
end
figure('Name', ['Dataset ', num2str(data_set)])
bar(mean_lost, 'k')
hold on;
errorbar(mean_lost, sem_lost, 'r+');
hold off;
set(gca, 'XTickLabel', bin_sizes)
title('percentage of spikes lost through binning')
xlabel('bin size (ms)')
ylabel('spikes lost (%)')
legend({'average over all cells', 'standard error of the mean'}, 'Location', 'NorthWest');