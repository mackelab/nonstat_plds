classdef preprocessor_alex < handle
    %DATA_ALEX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
%         path = '/Users/florian/ownCloud/Data';
        path = '/Users/mijung/code_mijung/Figure4/alex';
    end
    
    properties
        key;
        
        avg_corr;
        avg_instability;
        temporal_frequency;
        
        stim_duration;
        directions;
        
        stimulus_onset;
        conditions;
        spikes;
    end
    
    methods
        function self = preprocessor_alex(index)
            addpath(genpath(self.path));
            temp = load('ecker2014_subset.mat');
            
            % Save key, who knows if we ever need it.
            self.key = temp.keys(index);
            
            % Read data structure
            temp = temp.data(index);
            
            self.avg_corr = temp.avg_corr;
            self.avg_instability = temp.avg_instability;
            self.temporal_frequency = temp.temporal_frequency;
            
            self.stim_duration = temp.stim_duration;
            self.directions = pi * temp.directions / 180;
            
            self.stimulus_onset = temp.stimulus_onset;
            self.conditions = temp.conditions;
            self.spikes = temp.spikes;
            
        end
        
        function [x, y] = get_all_data(self, bin_size, binarize)
            % Inputs:
            %   bin_size: size of delta t in ms.
            %   binarize: Enable or disable binarization. (optional, default: true)
            %
            % Returns:
            %   x: (DxTxN)
            %   y: (CxTxN)
            
            if nargin < 3
                binarize = true;
            end
            
            N = size(self.stimulus_onset, 1);
            
            T = ceil((2 * self.stim_duration) / bin_size);
            
            onset = round(1000 / bin_size) + 1;
            offset = onset + round(self.stim_duration / bin_size);
            
            x = zeros(2, T, N);
            
            for trial = 1:N
                theta = self.directions(self.conditions(trial));
                x(1, onset:offset, trial) = cos(theta);
                x(2, onset:offset, trial) = sin(theta);
            end
            
            y = self.get_binned_spikes(bin_size, binarize);
        end
        
        function [x, y] = get_data_for_ori(self, orientation, bin_size, binarize)
            % Inputs:
            %   orientation: Orientation to use in degree
            %   bin_size: size of delta t in ms.
            %   binarize: Enable or disable binarization. (optional, default: true)
            %
            % Returns:
            %   x: (DxTxN)
            %   y: (CxTxN)
            
            if nargin < 4
                binarize = true;
            end
            
            % Generate x
            cond_index = find(self.directions == pi * orientation / 180);
            
            
            
            if length(cond_index) ~= 1
                error('Unknown or unprecise stimulus orientation supplied.')
            end
            trial_indices = self.conditions == cond_index;
            
            N = sum(trial_indices);
            
            T = ceil((2 * self.stim_duration) / bin_size);
            
            onset = round(1000 / bin_size) + 1;
            offset = onset - 1 + round(self.stim_duration / bin_size);
            
            x = zeros(3, T);
            
            max_phase = self.temporal_frequency * 2 * pi * (self.stim_duration / 1000);
            
            phase = 0:2 * max_phase / T :max_phase - 2 * max_phase / T;
            
            x(1, onset:offset) = 1;
            x(2, onset:offset) = cos(phase);
            x(3, onset:offset) = sin(phase);
            
            x = repmat(x, [1, 1, N]);
            
            % Generate y
            y = self.get_binned_spikes(bin_size, binarize, @(x) x(trial_indices, :));
        end
        
        function y = get_binned_spikes(self, bin_size, binarize, restrict_func)
            %GET_BINNED_SPIKES returns binned spikes.
            %   bin_size: Binning size in ms
            %   binarize: Enable or disable binarization.
            %   restrict_func: Function that allows you to restrict cell or trials to return. (Optional)
            %
            %   Returns a (CxTxN) matrix.
            %
            % IMPORTANT: If you want unbinarized spikes, just set binarize to false.
            
            if nargin < 4
                restrict_func = @(x) x;
            end
            
            restricted_spikes = restrict_func(self.spikes);
            
            [N, C] = size(restricted_spikes);
            
            % Trial length is shorted, to removes the last second of
            % recording.
            T = ceil(2 * self.stim_duration / bin_size);
            
            if binarize
                y = false(C, T, N);
            else
                y = zeros(C, T, N);
            end
            
            for unit = 1:C
                for trial = 1:N
                    spike_times = ceil((restricted_spikes{trial, unit} + 1000) / bin_size);
                    
                    % Remove last spikes, that reach into next trial
                    spike_times = spike_times(spike_times <= T);
                    
                    if binarize
                        y(unit, spike_times, trial) = true;
                    else
                        if(~isempty(spike_times))
                            % Find unique bins and their count.
                            unique_values = find([true; diff(spike_times) ~= 0; true]);
                            spike_times   = spike_times(unique_values(1:end-1));
                            spike_count   = diff(unique_values);
                            
                            y(unit, spike_times, trial) = spike_count;
                        end
                    end % if binarize
                end % for trial
            end % for unit
        end
    end
    
    methods (Static)
        function n = available_sessions()
            addpath(genpath(self.path));
            temp = load('ecker2014_subset.mat');
            n = length(temp.data);
        end
    end
end

