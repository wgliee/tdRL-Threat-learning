function results = decode_high_low_freezing_cs_function(traceeach_zs, freezing_each, clustidx_each, varargin)
%% Decode high-freezing CS vs. low-freezing CS
% This function classifies neural population activity during high-freezing
% and low-freezing CS trials using a linear discriminant classifier. It
% also runs a random trial-ranking control and compares test accuracy.
%
% @author: Mingyang Wei, Shanghai Jiao Tong University, 2026
%
% Inputs:
%       - traceeach_zs   : {mouse}(neurons x frames), z-scored neural traces
%       - freezing_each : {mouse}(frames x 1), freezing-state data
%       - clustidx_each  : {mouse}(neurons x 1), neuronal subgroup labels
%
% Optional name-value inputs:
%       - 'TraceEachZsDs'   : pre-downsampled traces used when DsTimes = 10.
%                             Default: {}
%       - 'DsTimes'         : downsampling factor. Default: 1
%       - 'SelectSubgroup'  : neuronal subgroup to include. Default: 4
%       - 'SelectThreshold' : fraction denominator used to select low and
%                             high trials. Default: 4
%       - 'NumFold'         : number of cross-validation folds. Default: 2
%                             Default: 1
%       - 'IterNum'         : number of repeated iterations. Default: 100
%       - 'DiscrimType'     : fitcdiscr discriminant type.
%                             Default: 'pseudolinear'
%       - 'PlotResults'     : whether to plot results. Default: true
%
% Outputs:
%       - results.savevar              : test accuracy for real and random data
%       - results.lowall               : mean low-freezing CS responses
%       - results.highall              : mean high-freezing CS responses
%       - results.lotrace              : low-freezing CS traces for each mouse
%       - results.hitrace              : high-freezing CS traces for each mouse
%       - results.p_accuracy           : real vs. random accuracy t-test p value
%       - results.p_response           : low vs. high response t-test p value
%       - results.params               : parameters used in the analysis
%

p = inputParser;
addRequired(p, 'traceeach_zs');
addRequired(p, 'freezing_each');
addRequired(p, 'clustidx_each');
addParameter(p, 'TraceEachZsDs', {});
addParameter(p, 'DsTimes', 1);
addParameter(p, 'SelectSubgroup', 4);
addParameter(p, 'SelectThreshold', 4);
addParameter(p, 'NumFold', 2);
addParameter(p, 'IterNum', 100);
addParameter(p, 'DiscrimType', 'pseudolinear');
addParameter(p, 'PlotResults', true);
parse(p, traceeach_zs, freezing_each, clustidx_each, varargin{:});

params = p.Results;
n = numel(traceeach_zs);

savevar = [];
lowall = [];
highall = [];
hitrace = cell(1, n);
lotrace = cell(1, n);
mice_v_acc_all = [];
mice_t_acc_all = [];
mean_mice_v_acc = zeros(1, 2);
mean_mice_t_acc = zeros(1, 2);

% Mean fluorescence and CS timing.
ds_times = params.DsTimes;
tone_start = [[180 330 450 630 720] * 10, ...
    [180 260 350 410 490 570 630 680 730 780 870 930] * 10 + 7600] ./ ds_times;
tone_end = ([[180 330 450 630 720] * 10, ...
    [180 260 350 410 490 570 630 680 730 780 870 930] * 10 + 7600] + 300) ./ ds_times;

include_length = 300 ./ ds_times;

if params.PlotResults
    figure;
    hold on;
end

for CreateRandomControl = 0:1
    mice_v_acc = [];
    mice_t_acc = [];
    foldwiseTestAcc = zeros(n, params.NumFold, params.IterNum); 

    for micenum = 1:n
        fprintf('\nNow processing micenum: %d', micenum)

        if ds_times == 10
            micetrace = params.TraceEachZsDs{micenum};
        elseif ds_times == 1
            micetrace = traceeach_zs{micenum};
        else
            micetrace = downsample_trace(traceeach_zs{micenum}, ds_times);
        end

        micefreeze = freezing_each{micenum};
        micefreeze = downsample_trace(micefreeze', ds_times)';

        % Select neuronal subpopulation.
        curmiceidx = clustidx_each{micenum};
        micetrace = micetrace(curmiceidx == params.SelectSubgroup, :);

        tone_dataset = [];
        freeze_dataset = [];

        % Prepare CS and freezing datasets.
        for tones = 1:size(tone_start, 2)
            tone_dataset = cat(3, tone_dataset, ...
                micetrace(:, tone_start(tones):tone_start(tones)+include_length-1));
            freeze_dataset = [freeze_dataset, ...
                micefreeze(tone_start(tones):tone_end(tones)-1)]; 
        end

        neuron_num = size(tone_dataset, 1);

        % Rank trials by freezing.
        trail_freeze = sum(freeze_dataset, 1);
        [~, sortidx] = sort(trail_freeze, 'ascend');

        validationAccuracy = zeros(1, params.IterNum);
        testAccuracy = zeros(1, params.IterNum);

        for iters = 1:params.IterNum
            % Randomly rank trials for the control analysis.
            if CreateRandomControl
                sortidx = randperm(max(sortidx));
            end

            select_trails = floor(size(tone_start, 2) ./ params.SelectThreshold);
            low_dataset_idx = sortidx(1:select_trails);
            high_dataset_idx = sortidx(size(tone_start, 2)-select_trails+1:end);

            % Display included freezing trials.
            if iters == 1
                [trail_freeze(low_dataset_idx), trail_freeze(high_dataset_idx)] 
            end

            low_freeze_dataset = tone_dataset(:, :, low_dataset_idx);
            high_freeze_dataset = tone_dataset(:, :, high_dataset_idx);

            % Save low- and high-freezing traces from the real data.
            if ~CreateRandomControl
                if iters == 1
                    lofreeze = squeeze(mean(low_freeze_dataset, 1))';
                    hifreeze = squeeze(mean(high_freeze_dataset, 1))';
                    lowall = [lowall; mean(lofreeze, 1)]; 
                    highall = [highall; mean(hifreeze, 1)]; 

                    lotrace{micenum} = reshape(permute(low_freeze_dataset, ...
                        [1, 3, 2]), neuron_num, [])';
                    hitrace{micenum} = reshape(permute(high_freeze_dataset, ...
                        [1, 3, 2]), neuron_num, [])';
                end
            end

            % Cross-validation partition by trial.
            cvlo = cvpartition(select_trails, "KFold", params.NumFold);
            cvhi = cvpartition(select_trails, "KFold", params.NumFold);

            valAcc = zeros(1, params.NumFold);
            tstAcc = zeros(1, params.NumFold);

            for k = 1:params.NumFold
                trainlo = training(cvlo, k);
                testlo = test(cvlo, k);
                trainhi = training(cvhi, k);
                testhi = test(cvhi, k);

                assert(isempty(find(trainlo .* testlo, 1)), ...
                    'Error: intersection found in training and test sets')
                assert(isempty(find(trainhi .* testhi, 1)), ...
                    'Error: intersection found in training and test sets')

                % Labelling: 1 = low freezing, 0 = high freezing.
                training_data_lo = reshape(permute( ...
                    low_freeze_dataset(:, :, trainlo), [1, 3, 2]), neuron_num, [])';
                training_data_hi = reshape(permute( ...
                    high_freeze_dataset(:, :, trainhi), [1, 3, 2]), neuron_num, [])';
                training_data = [training_data_lo; training_data_hi];
                training_label = [ones(size(training_data_lo, 1), 1); ...
                    zeros(size(training_data_hi, 1), 1)];

                test_data_lo = reshape(permute( ...
                    low_freeze_dataset(:, :, testlo), [1, 3, 2]), neuron_num, [])';
                test_data_hi = reshape(permute( ...
                    high_freeze_dataset(:, :, testhi), [1, 3, 2]), neuron_num, [])';
                test_data = [test_data_lo; test_data_hi];
                test_label = [ones(size(test_data_lo, 1), 1); ...
                    zeros(size(test_data_hi, 1), 1)];

                model = fitcdiscr(training_data, training_label, ...
                    'DiscrimType', params.DiscrimType);

                % Validation accuracy.
                valPredictions = predict(model, training_data);
                valAcc(k) = sum(training_label == valPredictions) / length(training_label);

                % Test accuracy.
                testPredictions = predict(model, test_data);
                tstAcc(k) = sum(test_label == testPredictions) / length(test_label);
            end

            validationAccuracy(iters) = mean(valAcc, "all", "omitnan");
            testAccuracy(iters) = mean(tstAcc, "all", "omitnan");
        end

        mice_v_acc = vertcat(mice_v_acc, mean(validationAccuracy));
        mice_t_acc = vertcat(mice_t_acc, mean(testAccuracy));
    end

    mean_mice_v_acc(CreateRandomControl + 1) = mean(mice_v_acc);
    mean_mice_t_acc(CreateRandomControl + 1) = mean(mice_t_acc);
    mice_v_acc_all(:, CreateRandomControl + 1) = mice_v_acc;
    mice_t_acc_all(:, CreateRandomControl + 1) = mice_t_acc;

    if params.PlotResults
        rowdata_columnumber_plot(mice_t_acc);
    end

    savevar(:, :, CreateRandomControl + 1) = mice_t_acc;
end

[h_accuracy, p_accuracy] = ttest(mean(savevar(:, :, 1), 2), ...
    mean(savevar(:, :, 2), 2));

if params.PlotResults
    ylabel('Average Test Accuracy');
    title(strcat('Test Accuracy, p=', string(p_accuracy)));
    legend("Real Data", "Random Control");
    hold off;
end

% Inspect high- and low-freezing CS responses.
[~, p_response] = ttest(sum(lowall, 2), sum(highall, 2));

if params.PlotResults
    figure;
    hold on;
    rowdata_columnumber_plot(lowall);
    rowdata_columnumber_plot(highall);
    hold off;
    legend("low freezing", "high freezing");
    title(string(p_response));
end


params_str = sprintf([ ...
    'Mean Fluorescence\n' ...
    'ds_times = %d\n' ...
    'Tone Start: %s\n' ...
    'Tone End: %s\n' ...
    'Include Length = %d\n' ...
    'numFold = %d\n' ...
    'iternum = %d\n'], ...
    ds_times, mat2str(tone_start), mat2str(tone_end), ...
    include_length, params.NumFold, params.IterNum);

results.savevar = savevar;
results.lowall = lowall;
results.highall = highall;
results.lotrace = lotrace;
results.hitrace = hitrace;
results.mice_v_acc = mice_v_acc_all;
results.mice_t_acc = mice_t_acc_all;
results.mean_mice_v_acc = mean_mice_v_acc;
results.mean_mice_t_acc = mean_mice_t_acc;
results.h_accuracy = h_accuracy;
results.p_accuracy = p_accuracy;
results.p_response = p_response;
results.tone_start = tone_start;
results.tone_end = tone_end;
results.include_length = include_length;
results.params_str = params_str;
results.params = params;

end
