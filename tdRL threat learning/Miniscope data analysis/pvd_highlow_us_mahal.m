function results = pvd_highlow_us_mahal(traceeach_zs, freezing_each, clustidx_each, varargin)
%% PVD distance from high/low-freezing CS to all US
% This function calculates normalized Mahalanobis distances from population
% activity during low- and high-freezing CS trials to population activity
% during all US trials. A random trial-order control is also calculated.
%
% @author: Mingyang Wei, Shanghai Jiao Tong University, 2026
%
% Inputs:
%       - traceeach_zs   : {mouse}(neurons x frames), z-scored neural traces
%       - freezing_each : {mouse}, frame-by-frame freezing data
%       - clustidx_each  : {mouse}, neuronal subgroup index
%
% Optional name-value inputs:
%       - 'DsTimes'            : downsampling factor. Default: 1
%       - 'SelectSubgroup'     : neuronal subgroup. Default: 4
%       - 'SelectThreshold'    : denominator for selecting high/low trials.
%                                Default: 4
%       - 'ShockIncludeLength' : US window length before downsampling.
%                                Default: 40
%       - 'IntervalLength'     : freezing averaging window before
%                                downsampling. Default: 150
%       - 'IncludeLength'      : CS window length before downsampling.
%                                Default: 300
%       - 'ControlIterNum'     : random-control iterations. Default: 100
%       - 'PlotResults'        : whether to plot results. Default: true
%

p = inputParser;
addRequired(p, 'traceeach_zs');
addRequired(p, 'freezing_each');
addRequired(p, 'clustidx_each');
addParameter(p, 'DsTimes', 1);
addParameter(p, 'SelectSubgroup', 4);
addParameter(p, 'SelectThreshold', 4);
addParameter(p, 'ShockIncludeLength', 40);
addParameter(p, 'IntervalLength', 150);
addParameter(p, 'IncludeLength', 300);
addParameter(p, 'ControlIterNum', 100);
addParameter(p, 'PlotResults', true);
parse(p, traceeach_zs, freezing_each, clustidx_each, varargin{:});

params = p.Results;
n = numel(traceeach_zs);

ds_times = params.DsTimes;
select_subgroup = params.SelectSubgroup;
select_threshold = params.SelectThreshold;
shock_include_length = params.ShockIncludeLength ./ ds_times;
interval_length = params.IntervalLength ./ ds_times;
include_length = params.IncludeLength ./ ds_times;
control_iter_num = params.ControlIterNum;

shock_start = ([180 330 450 630 720] * 10 + 300) ./ ds_times;
tone_start = ([[180 330 450 630 720] * 10, ...
    [180 260 350 410 490 570 630 680 730 780 870 930] * 10 + 7600]) ./ ds_times;
tone_end = ([[180 330 450 630 720] * 10, ...
    [180 260 350 410 490 570 630 680 730 780 870 930] * 10 + 7600] + 300) ./ ds_times;

pvddist_low2early = [];
pvddist_high2early = [];
pvddist_low2early_ctrl = [];
pvddist_high2early_ctrl = [];
hi_low_trials_used = [];

for CreateRdmCtrl = 0:1
    if CreateRdmCtrl == 0
        for micenum = 1:n
            fprintf('\nNow processing micenum: %d', micenum)

            if ds_times == 1
                micetrace = traceeach_zs{micenum};
            else
                micetrace = downsample_trace(traceeach_zs{micenum}, ds_times);
            end

            micefreeze = freezing_each{micenum};
            micefreeze = downsample_trace(micefreeze', ds_times)';
            freezing_probability = timeSeriesAverage(micefreeze, interval_length)';

            curmiceidx = clustidx_each{micenum};
            micetrace = micetrace(curmiceidx == select_subgroup, :);

            tone_dataset = [];
            freeze_dataset = [];
            shock_dataset = [];
            freezeprob_dataset = [];

            for shocks = 1:size(shock_start, 2)
                shock_dataset = cat(3, shock_dataset, ...
                    micetrace(:, shock_start(shocks):shock_start(shocks)+shock_include_length-1));
            end

            for tones = 1:size(tone_start, 2)
                tone_dataset = cat(3, tone_dataset, ...
                    micetrace(:, tone_start(tones):tone_start(tones)+include_length-1));
                freeze_dataset = [freeze_dataset, ...
                    micefreeze(tone_start(tones):tone_end(tones)-1)];
                freezeprob_dataset = [freezeprob_dataset, ...
                    freezing_probability(tone_start(tones):tone_end(tones)-1)];
            end

            trail_freeze = sum(freeze_dataset, 1);
            [~, sortidx] = sort(trail_freeze, 'ascend');

            select_trails = floor(size(tone_start, 2) ./ select_threshold);
            low_dataset_idx = sortidx(1:select_trails);
            high_dataset_idx = sortidx(size(tone_start, 2)-select_trails+1:end);
            hi_low_trials_used(:, :, micenum) = [high_dataset_idx; low_dataset_idx];

            low_freeze_dataset = tone_dataset(:, :, low_dataset_idx);
            high_freeze_dataset = tone_dataset(:, :, high_dataset_idx);

            low_freezing_prob = reshape(freezeprob_dataset(:, low_dataset_idx), [], 1); 
            high_freezing_prob = reshape(freezeprob_dataset(:, high_dataset_idx), [], 1);

            lo_freeze_mahal = reshape(permute(low_freeze_dataset, [3 2 1]), ...
                [], size(low_freeze_dataset, 1));
            high_freeze_mahal = reshape(permute(high_freeze_dataset, [3 2 1]), ...
                [], size(high_freeze_dataset, 1));

            mahal_ref_US = shock_dataset(:, :, 1:5);
            reshaped_early = permute(mahal_ref_US, [3 2 1]);
            mahal_ref_US = reshape(mahal_ref_US, [], size(mahal_ref_US, 1));

            pvddist_low2early(micenum, 1) = ...
                mean(mahal(lo_freeze_mahal, mahal_ref_US), 'all') ./ ...
                mean(mahal(mahal_ref_US, mahal_ref_US), 'all');
            pvddist_high2early(micenum, 1) = ...
                mean(mahal(high_freeze_mahal, mahal_ref_US), 'all') ./ ...
                mean(mahal(mahal_ref_US, mahal_ref_US), 'all');
        end

    elseif CreateRdmCtrl == 1
        for micenum = 1:n
            fprintf('\nNow processing micenum: %d', micenum)

            if ds_times == 1
                micetrace = traceeach_zs{micenum};
            else
                micetrace = downsample_trace(traceeach_zs{micenum}, ds_times);
            end

            micefreeze = freezing_each{micenum};
            micefreeze = downsample_trace(micefreeze', ds_times)';
            freezing_probability = timeSeriesAverage(micefreeze, interval_length)';

            curmiceidx = clustidx_each{micenum};
            micetrace = micetrace(curmiceidx == select_subgroup, :);

            tone_dataset = [];
            freeze_dataset = [];
            shock_dataset = [];
            freezeprob_dataset = [];

            for shocks = 1:size(shock_start, 2)
                shock_dataset = cat(3, shock_dataset, ...
                    micetrace(:, shock_start(shocks):shock_start(shocks)+shock_include_length-1));
            end

            for tones = 1:size(tone_start, 2)
                tone_dataset = cat(3, tone_dataset, ...
                    micetrace(:, tone_start(tones):tone_start(tones)+include_length-1));
                freeze_dataset = [freeze_dataset, ...
                    micefreeze(tone_start(tones):tone_end(tones)-1)];
                freezeprob_dataset = [freezeprob_dataset, ...
                    freezing_probability(tone_start(tones):tone_end(tones)-1)];
            end

            trail_freeze = sum(freeze_dataset, 1);
            [~, sortidx] = sort(trail_freeze, 'ascend');

            for iters = 1:control_iter_num
                sortidx = sortidx(randperm(length(sortidx)));
                select_trails = floor(size(tone_start, 2) ./ select_threshold);
                low_dataset_idx = sortidx(1:select_trails);
                high_dataset_idx = sortidx(size(tone_start, 2)-select_trails+1:end);

                low_freeze_dataset = tone_dataset(:, :, low_dataset_idx);
                high_freeze_dataset = tone_dataset(:, :, high_dataset_idx);

                low_freezing_prob = reshape(freezeprob_dataset(:, low_dataset_idx), [], 1); 
                high_freezing_prob = reshape(freezeprob_dataset(:, high_dataset_idx), [], 1); 

                lo_freeze_mahal = reshape(permute(low_freeze_dataset, [3 2 1]), ...
                    [], size(low_freeze_dataset, 1));
                high_freeze_mahal = reshape(permute(high_freeze_dataset, [3 2 1]), ...
                    [], size(high_freeze_dataset, 1));

                mahal_ref_US = shock_dataset(:, :, 1:5);
                reshaped_early = permute(mahal_ref_US, [3 2 1]); 
                mahal_ref_US = reshape(mahal_ref_US, [], size(mahal_ref_US, 1));

                lowtemp(iters, 1) = ...
                    mean(mahal(lo_freeze_mahal, mahal_ref_US), 'all') ./ ...
                    mean(mahal(mahal_ref_US, mahal_ref_US), 'all');
                hightemp(iters, 1) = ...
                    mean(mahal(high_freeze_mahal, mahal_ref_US), 'all') ./ ...
                    mean(mahal(mahal_ref_US, mahal_ref_US), 'all');
            end

            pvddist_low2early_ctrl(micenum, 1) = mean(lowtemp, 'all');
            pvddist_high2early_ctrl(micenum, 1) = mean(hightemp, 'all');
        end
    end
end

pvd_results = [pvddist_low2early, pvddist_high2early, ...
    pvddist_low2early_ctrl, pvddist_high2early_ctrl];
variableNames = {'Low2US', 'High2US', 'Ctrl_Low2US', 'Ctrl_High2US'};

p_values = zeros(4, 4);
for i = 1:3
    for j = i+1:4
        [~, p_values(i, j)] = ttest(pvd_results(:, i), pvd_results(:, j));
    end
end

if params.PlotResults
    rowdata_columnumber_plot(pvd_results)
    xticks(1:4)
    xticklabels(variableNames)
end

results.pvd_results = pvd_results;
results.p_values = p_values;
results.hi_low_trials_used = hi_low_trials_used;
results.pvddist_low2early = pvddist_low2early;
results.pvddist_high2early = pvddist_high2early;
results.pvddist_low2early_ctrl = pvddist_low2early_ctrl;
results.pvddist_high2early_ctrl = pvddist_high2early_ctrl;
results.variableNames = variableNames;
results.params = params;

end
