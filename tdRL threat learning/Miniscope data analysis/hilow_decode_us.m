function results = hilow_decode_us( ...
    traceeach_zs, traceeach_zs_ds, traceeach_S, freezing_each, clustidx_each, varargin)
%% Decode shock as high-freezing CS with matched subgroup neuron numbers
% This function compares subgroup 1 and subgroup 4 after randomly matching
% their neuron numbers. For each US trial, the corresponding CS trial is
% excluded before selecting low- and high-freezing CS trials. A classifier
% is trained on low- versus high-freezing CS activity and used to calculate
% the proportion of shock frames predicted as high freezing.
%
% Inputs:
%       - traceeach_zs:    {mouse}(neurons x frames), z-scored traces
%       - traceeach_zs_ds: {mouse}(neurons x frames), downsampled traces
%       - traceeach_S:     {mouse}(neurons x frames), neural traces used when
%                          DsTimes = 1, following the original script
%       - freezing_each:   {mouse}, frame-wise freezing data
%       - clustidx_each:   {mouse}, subgroup label for each neuron
%
% Optional name-value inputs:
%       - 'DsTimes'         : downsampling factor. Default: 1
%       - 'SelectSubgroups' : subgroup labels to compare. Default: [1 4]
%       - 'SelectThreshold' : fraction used to select low/high CS trials.
%                             Default: 4
%       - 'RealIterNum'     : real decoding iterations. Default: 100
%       - 'ChanceIterNum'   : random-control iterations. Default: 100
%       - 'DiscrimType'     : fitcdiscr discriminant type.
%                             Default: 'pseudolinear'
%       - 'PlotResults'     : whether to plot final mouse-level results.
%                             Default: true
%
% Outputs:
%       - results.pred_shock_as_high_real
%       - results.pred_shock_as_high_chance
%       - results.real_mouse_mean_pred
%       - results.chance_mouse_mean_pred
%       - results.real_mean_pred_by_subgroup
%       - results.chance_mean_pred_by_subgroup
%       - results.real_minus_chance_by_subgroup
%       - results.p_real_vs_chance
%       - results.p_subgroup_real
%       - results.p_subgroup_chance
%
% Required helper functions:
%       - random_subsample_getidx
%       - downsample_trace
%       - rowdata_columnumber_plot

p = inputParser;
addRequired(p, 'traceeach_zs');
addRequired(p, 'traceeach_zs_ds');
addRequired(p, 'traceeach_S');
addRequired(p, 'freezing_each');
addRequired(p, 'clustidx_each');
addParameter(p, 'DsTimes', 1);
addParameter(p, 'SelectSubgroups', [1 4]);
addParameter(p, 'SelectThreshold', 4);
addParameter(p, 'RealIterNum', 100);
addParameter(p, 'ChanceIterNum', 100);
addParameter(p, 'DiscrimType', 'pseudolinear');
addParameter(p, 'PlotResults', true);
parse(p, traceeach_zs, traceeach_zs_ds, traceeach_S, ...
    freezing_each, clustidx_each, varargin{:});

params = p.Results;
n = numel(freezing_each);

ds_times = params.DsTimes;
select_subgroups = params.SelectSubgroups;
nSubgroups = numel(select_subgroups);
select_threshold = params.SelectThreshold;

shock_include_length = 100 ./ ds_times;

tone_start = ([[180 330 450 630 720] * 10, ...
    [180 260 350 410 490 570 630 680 730 780 870 930] * 10 + 7600]) ./ ds_times;
tone_end = ([[180 330 450 630 720] * 10, ...
    [180 260 350 410 490 570 630 680 730 780 870 930] * 10 + 7600] + 300) ./ ds_times;

shock_start = ([[180 330 450 630 720] * 10] + 300) ./ ds_times;
shock_end = ([[180 330 450 630 720] * 10] + 300 + shock_include_length) ./ ds_times;

include_length = 300 ./ ds_times;
nCS = numel(tone_start);
nUS = numel(shock_start);
all_cs_idx = 1:nCS;
all_us_idx = 1:nUS;

for CreateRandomControl = 0:1
    if CreateRandomControl
        iternum = params.ChanceIterNum;
    else
        iternum = params.RealIterNum;
    end

    % Dimensions: mouse x subgroup x US x iteration.
    pred_shock_as_high = nan(n, nSubgroups, nUS, iternum);

    for micenum = 1:n
        fprintf('\nNow processing micenum: %d', micenum)

        if ds_times == 10
            micetrace_all = traceeach_S{micenum};
        elseif ds_times == 1
            micetrace_all = traceeach_S{micenum};
        else
            micetrace_all = downsample_trace(traceeach_S{micenum}, ds_times);
        end

        micefreeze = freezing_each{micenum};
        micefreeze = downsample_trace(micefreeze', ds_times)';

        curmiceidx = clustidx_each{micenum};
        micetrace_subgroup_full = cell(1, nSubgroups);

        for subgroup_i = 1:nSubgroups
            subgroup_label = select_subgroups(subgroup_i);
            micetrace_subgroup_full{subgroup_i} = ...
                micetrace_all(curmiceidx == subgroup_label, :);
        end

        % Rank CS trials by freezing.
        freeze_dataset = [];
        for tones = 1:nCS
            freeze_dataset = [freeze_dataset, ...
                micefreeze(tone_start(tones):tone_end(tones) - 1)];
        end

        trail_freeze = sum(freeze_dataset, 1);
        [~, sortidx_by_freezing_all] = sort(trail_freeze, 'ascend');

        % Prepare subgroup-specific CS and US datasets.
        tone_dataset_full = cell(1, nSubgroups);
        shock_dataset_full = cell(1, nSubgroups);

        for subgroup_i = 1:nSubgroups
            micetrace_tmp = micetrace_subgroup_full{subgroup_i};

            tone_dataset_tmp = [];
            shock_dataset_tmp = [];

            for tones = 1:nCS
                tone_dataset_tmp = cat(3, tone_dataset_tmp, ...
                    micetrace_tmp(:, tone_start(tones):tone_start(tones) + include_length - 1));
            end

            for shocks = 1:nUS
                shock_dataset_tmp = cat(3, shock_dataset_tmp, ...
                    micetrace_tmp(:, shock_start(shocks):shock_end(shocks) - 1));
            end

            tone_dataset_full{subgroup_i} = tone_dataset_tmp;
            shock_dataset_full{subgroup_i} = shock_dataset_tmp;
        end

        for iters = 1:iternum
            % Match subgroup 1 and subgroup 4 neuron numbers.
            [idx_subgroup1, idx_subgroup4] = random_subsample_getidx( ...
                micetrace_subgroup_full{1}, micetrace_subgroup_full{2});
            subsample_idx_this_iter = {idx_subgroup1, idx_subgroup4};

            for subgroup_i = 1:nSubgroups
                neuron_keep_idx = subsample_idx_this_iter{subgroup_i};

                tone_dataset = tone_dataset_full{subgroup_i}(neuron_keep_idx, :, :);
                shock_dataset = shock_dataset_full{subgroup_i}(neuron_keep_idx, :, :);
                neuron_num = size(tone_dataset, 1);

                for test_us_id = all_us_idx
                    % Leave the matched CS trial out before selecting low/high CS.
                    candidate_cs_idx = setdiff(all_cs_idx, test_us_id, 'stable');

                    if CreateRandomControl
                        sortidx_current = candidate_cs_idx(randperm(numel(candidate_cs_idx)));
                    else
                        sortidx_current = sortidx_by_freezing_all( ...
                            ~ismember(sortidx_by_freezing_all, test_us_id));
                    end

                    select_trails = floor(numel(candidate_cs_idx) ./ select_threshold);
                    low_dataset_idx = sortidx_current(1:select_trails);
                    high_dataset_idx = sortidx_current( ...
                        numel(sortidx_current) - select_trails + 1:end);

                    low_freeze_dataset = tone_dataset(:, :, low_dataset_idx);
                    high_freeze_dataset = tone_dataset(:, :, high_dataset_idx);

                    % Train classifier: 1 = low freezing; 0 = high freezing.
                    training_data_lo = reshape(permute( ...
                        low_freeze_dataset, [1 3 2]), neuron_num, [])';
                    training_data_hi = reshape(permute( ...
                        high_freeze_dataset, [1 3 2]), neuron_num, [])';
                    training_data = [training_data_lo; training_data_hi];
                    training_label = [ones(size(training_data_lo, 1), 1); ...
                        zeros(size(training_data_hi, 1), 1)];

                    model = fitcdiscr(training_data, training_label, ...
                        'DiscrimType', params.DiscrimType);

                    % Predict the current shock trial as high freezing.
                    shock_dataset_for_test = shock_dataset(:, :, test_us_id);
                    test_data_shock = reshape(permute( ...
                        shock_dataset_for_test, [1 3 2]), neuron_num, [])';

                    pred_shock_as_high(micenum, subgroup_i, test_us_id, iters) = ...
                        sum(predict(model, test_data_shock) == 0) ./ ...
                        size(test_data_shock, 1);
                end
            end
        end
    end

    if CreateRandomControl
        pred_shock_as_high_chance = pred_shock_as_high;
    else
        pred_shock_as_high_real = pred_shock_as_high;
    end
end

% Average across iterations and then across US trials.
pred_shock_as_high_real_mean = mean( ...
    pred_shock_as_high_real, 4, 'omitnan');
pred_shock_as_high_chance_mean = mean( ...
    pred_shock_as_high_chance, 4, 'omitnan');

real_mouse_mean_pred = squeeze(mean( ...
    pred_shock_as_high_real_mean, 3, 'omitnan'));
chance_mouse_mean_pred = squeeze(mean( ...
    pred_shock_as_high_chance_mean, 3, 'omitnan'));

% Final subgroup results.
real_mean_pred_by_subgroup = mean(real_mouse_mean_pred, 1, 'omitnan');
chance_mean_pred_by_subgroup = mean(chance_mouse_mean_pred, 1, 'omitnan');
real_minus_chance_by_subgroup = ...
    real_mean_pred_by_subgroup - chance_mean_pred_by_subgroup;

% Paired real-vs-chance tests within each subgroup.
h_real_vs_chance = nan(1, nSubgroups);
p_real_vs_chance = nan(1, nSubgroups);

for subgroup_i = 1:nSubgroups
    [h_real_vs_chance(subgroup_i), p_real_vs_chance(subgroup_i)] = ...
        ttest(real_mouse_mean_pred(:, subgroup_i), ...
        chance_mouse_mean_pred(:, subgroup_i));
end

% Paired subgroup comparison after neuron-number matching.
[h_subgroup_real, p_subgroup_real] = ...
    ttest(real_mouse_mean_pred(:, 1), real_mouse_mean_pred(:, 2));
[h_subgroup_chance, p_subgroup_chance] = ...
    ttest(chance_mouse_mean_pred(:, 1), chance_mouse_mean_pred(:, 2));

if params.PlotResults
    figure;
    rowdata_columnumber_plot([ ...
        real_mouse_mean_pred(:, 1), chance_mouse_mean_pred(:, 1), ...
        real_mouse_mean_pred(:, 2), chance_mouse_mean_pred(:, 2)]);
    title(sprintf('Real vs chance, subgroup %d and subgroup %d', ...
        select_subgroups(1), select_subgroups(2)));
    xticklabels({ ...
        sprintf('SG%d real', select_subgroups(1)), ...
        sprintf('SG%d chance', select_subgroups(1)), ...
        sprintf('SG%d real', select_subgroups(2)), ...
        sprintf('SG%d chance', select_subgroups(2))});
end

fprintf('\n\nSubgroup labels: %s', mat2str(select_subgroups));
fprintf('\nReal mean pred by subgroup: %s', ...
    mat2str(real_mean_pred_by_subgroup, 4));
fprintf('\nChance mean pred by subgroup: %s', ...
    mat2str(chance_mean_pred_by_subgroup, 4));
fprintf('\nReal-vs-chance p by subgroup: %s', ...
    mat2str(p_real_vs_chance, 4));
fprintf('\nSubgroup real comparison p: %.4g', p_subgroup_real);
fprintf('\nSubgroup chance comparison p: %.4g\n', p_subgroup_chance);

results.pred_shock_as_high_real = pred_shock_as_high_real;
results.pred_shock_as_high_chance = pred_shock_as_high_chance;
results.pred_shock_as_high_real_mean = pred_shock_as_high_real_mean;
results.pred_shock_as_high_chance_mean = pred_shock_as_high_chance_mean;
results.real_mouse_mean_pred = real_mouse_mean_pred;
results.chance_mouse_mean_pred = chance_mouse_mean_pred;
results.real_mean_pred_by_subgroup = real_mean_pred_by_subgroup;
results.chance_mean_pred_by_subgroup = chance_mean_pred_by_subgroup;
results.real_minus_chance_by_subgroup = real_minus_chance_by_subgroup;
results.h_real_vs_chance = h_real_vs_chance;
results.p_real_vs_chance = p_real_vs_chance;
results.h_subgroup_real = h_subgroup_real;
results.p_subgroup_real = p_subgroup_real;
results.h_subgroup_chance = h_subgroup_chance;
results.p_subgroup_chance = p_subgroup_chance;
results.select_subgroups = select_subgroups;
results.tone_start = tone_start;
results.shock_start = shock_start;
results.params = params;

end
