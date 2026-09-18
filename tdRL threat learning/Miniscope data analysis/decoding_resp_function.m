function results = decoding_resp_function(traceeach_zs, decode_event, varargin)
%% Decode event vs. baseline response
% This function decodes neural population activity during tone or shock
%       periods from the corresponding pre-event baseline using a linear
%       discriminant classifier. It also generates a random-label control
%       and calculates accuracy, F1 score, and AUC.
%
% @author:  Mingyang Wei, Shanghai Jiao Tong University, 2026
%
% Inputs:
%       - traceeach_zs: {mouse}(neurons x frames), z-scored neural traces
%       - decode_event: event to decode from baseline
%             "tone"  : tone onset vs. pre-tone baseline
%             "shock" : shock onset vs. pre-shock baseline
%             "both"  : run both tone and shock decoding
%
% Optional name-value inputs:
%       - 'DsTimes'              : downsampling factor. Default: 1
%       - 'NumFold'              : number of CV folds. Default: 2
%       - 'IterNum'              : number of repeated CV iterations. Default: 100
%       - 'PlotResults'          : whether to plot test accuracy. Default: true
%       - 'CreateRandomControl'  : controls to run. Default: [0 1]
%       - 'ToneStartFrames'      : tone onset frames before downsampling.
%                                  Default: [180 330 450 630 720] * 10
%       - 'ShockStartFrames'     : shock onset frames before downsampling.
%                                  Default: [180 330 450 630 720] * 10 + 300
%       - 'ToneIncludeLength'    : tone analysis window before downsampling.
%                                  Default: 100
%       - 'ShockIncludeLength'   : shock analysis window before downsampling.
%                                  Default: 40
%       - 'DiscrimType'          : fitcdiscr discriminant type.
%                                  Default: 'pseudolinear'
%


if nargin < 2 || isempty(decode_event)
    decode_event = "tone";
end

decode_event = lower(string(decode_event));

p = inputParser;
addRequired(p, 'traceeach_zs');
addRequired(p, 'decode_event');
addParameter(p, 'DsTimes', 1);
addParameter(p, 'NumFold', 2);
addParameter(p, 'IterNum', 100);
addParameter(p, 'PlotResults', true);
addParameter(p, 'CreateRandomControl', [0 1]);
addParameter(p, 'ToneStartFrames', [180 330 450 630 720] * 10);
addParameter(p, 'ShockStartFrames', [180 330 450 630 720] * 10 + 300);
addParameter(p, 'ToneIncludeLength', 100);
addParameter(p, 'ShockIncludeLength', 40);
addParameter(p, 'DiscrimType', 'pseudolinear');
parse(p, traceeach_zs, decode_event, varargin{:});

params = p.Results;
n = numel(traceeach_zs);

if decode_event == "both"
    results.tone = run_single_event_decode("tone", traceeach_zs, n, params);
    results.shock = run_single_event_decode("shock", traceeach_zs, n, params);
elseif decode_event == "tone" || decode_event == "shock"
    results = run_single_event_decode(decode_event, traceeach_zs, n, params);
else
    error('decode_event must be "tone", "shock", or "both".');
end

end

function results = run_single_event_decode(decode_event, traceeach_zs, n, params)
%% Select event timing
if decode_event == "tone"
    event_start = params.ToneStartFrames ./ params.DsTimes;
    include_length = params.ToneIncludeLength ./ params.DsTimes;
elseif decode_event == "shock"
    event_start = params.ShockStartFrames ./ params.DsTimes;
    include_length = params.ShockIncludeLength ./ params.DsTimes;
else
    error('decode_event must be "tone" or "shock".');
end

% Keep frame indices as integers.
event_start = round(event_start);
include_length = round(include_length);

savevar = [];
f1Scores = [];
aucScores = [];
foldwiseTestAcc_all = [];
mice_v_acc_all = [];
mice_t_acc_all = [];

if params.PlotResults
    figure;
    hold on;
end

for CreateRandomControl = params.CreateRandomControl
    mice_v_acc = [];
    mice_t_acc = [];
    foldwiseTestAcc = zeros(n, numel(event_start), params.IterNum);
    foldwiseF1Score = nan(n, params.IterNum);
    foldwiseAUC = nan(n, params.IterNum);
    
    for micenum = 1:n
        if params.DsTimes == 1
            micetrace = traceeach_zs{micenum};
        else
            micetrace = downsample_trace(traceeach_zs{micenum}, params.DsTimes);
        end
        
        event_dataset = [];
        pre_event_dataset = [];
        
        % Prepare event and pre-event baseline datasets.
        for trial = 1:numel(event_start)
            event_dataset = cat(3, event_dataset, ...
                micetrace(:, event_start(trial):event_start(trial)+include_length-1));
            pre_event_dataset = cat(3, pre_event_dataset, ...
                micetrace(:, event_start(trial)-include_length:event_start(trial)-1));
        end
        
        neuron_num = size(event_dataset, 1);
        
        % Start training.
        validationAccuracy = zeros(1, params.IterNum);
        testAccuracy = zeros(1, params.IterNum);
        
        for iters = 1:params.IterNum
            cv = cvpartition(numel(event_start), "KFold", params.NumFold);
            valAcc = zeros(1, params.NumFold);
            tstAcc = zeros(1, params.NumFold);
            
            for k = 1:params.NumFold
                trainIdx = training(cv, k);
                testIdx = test(cv, k);
                
                % Labelling: 1 = event, 0 = pre-event baseline.
                training_data_event = reshape(permute(event_dataset(:,:,trainIdx), [1, 3, 2]), neuron_num, [])';
                training_data_pre_event = reshape(permute(pre_event_dataset(:,:,trainIdx), [1, 3, 2]), neuron_num, [])';
                training_data = [training_data_event; training_data_pre_event];
                training_label = [ones(size(training_data_event, 1), 1); zeros(size(training_data_pre_event, 1), 1)];
                
                test_data_event = reshape(permute(event_dataset(:,:,testIdx), [1, 3, 2]), neuron_num, [])';
                test_data_pre_event = reshape(permute(pre_event_dataset(:,:,testIdx), [1, 3, 2]), neuron_num, [])';
                test_data = [test_data_event; test_data_pre_event];
                test_label = [ones(size(test_data_event, 1), 1); zeros(size(test_data_pre_event, 1), 1)];
                
                if CreateRandomControl
                    training_label = training_label(randperm(length(training_label)));
                end
                
                model = fitcdiscr(training_data, training_label, 'DiscrimType', params.DiscrimType);
                
                % Validation accuracy.
                valPredictions = predict(model, training_data);
                valAcc(k) = sum(training_label == valPredictions) / length(training_label);
                
                % Test accuracy.
                [testPredictions, score] = predict(model, test_data);
                tstAcc(k) = sum(test_label == testPredictions) / length(test_label);
                foldwiseTestAcc(micenum, find(testIdx), iters) = tstAcc(k);
                
                % F1 score.
                [C, ~] = confusionmat(test_label, testPredictions);
                if size(C, 1) > 1
                    TP = C(2, 2);
                    FP = C(1, 2);
                    FN = C(2, 1);
                    Precision = TP / (TP + FP);
                    Recall = TP / (TP + FN);
                    F1 = 2 * ((Precision * Recall) / (Precision + Recall));
                else
                    F1 = NaN;
                end
                foldwiseF1Score(micenum, iters) = F1;
                
                % AUC.
                [~, ~, ~, AUC] = perfcurve(test_label, score(:,2), 1);
                foldwiseAUC(micenum, iters) = AUC;
            end
            
            % Store average accuracy for each iteration.
            validationAccuracy(iters) = mean(valAcc, "all", "omitnan");
            testAccuracy(iters) = mean(tstAcc, "all", "omitnan");
        end
        
        mice_v_acc = vertcat(mice_v_acc, mean(validationAccuracy));
        mice_t_acc = vertcat(mice_t_acc, mean(testAccuracy));
    end
    
    control_index = CreateRandomControl + 1;
    savevar(:,:,control_index) = mice_t_acc;
    f1Scores(:,:,control_index) = foldwiseF1Score;
    aucScores(:,:,control_index) = foldwiseAUC;
    foldwiseTestAcc_all(:,:,:,control_index) = foldwiseTestAcc;
    mice_v_acc_all(:,control_index) = mice_v_acc;
    mice_t_acc_all(:,control_index) = mice_t_acc;
    
    if params.PlotResults
        rowdata_columnumber_plot_omitnan(mice_t_acc);
    end
end

if params.PlotResults
    xlabel('Fold Number');
    ylabel('Average Test Accuracy');
    title(strcat('Test Accuracy Per ', upper(decode_event), '/Baseline as Test set'));
    legend("Real Data", "Random Control");
    xticks(1:numel(event_start));
    xticklabels(compose('%s %d', upper(decode_event), 1:numel(event_start)));
    hold off;
end

results.savevar = savevar;
results.f1Scores = f1Scores;
results.aucScores = aucScores;
results.foldwiseTestAcc = foldwiseTestAcc_all;
results.mice_v_acc = mice_v_acc_all;
results.mice_t_acc = mice_t_acc_all;
results.event_start = event_start;
results.include_length = include_length;
results.decode_event = decode_event;
results.params = params;

end
