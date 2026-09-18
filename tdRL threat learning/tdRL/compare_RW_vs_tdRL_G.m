%% Compare R-W and tdRL predictions at the five measured trials
% Run RW_G.m and tdRL_G.m first, or use run_example.m to execute the full
% workflow in the correct order.

scriptDir = fileparts(mfilename('fullpath'));
if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = fullfile(scriptDir, 'results');
end

rwFile = fullfile(outputDir, 'RW_fit_results.mat');
tdrlFile = fullfile(outputDir, 'tdRL_fit_results.mat');
if ~isfile(rwFile) || ~isfile(tdrlFile)
    error('ModelComparison:MissingFitResults', ...
        ['Fit-result files were not found in "%s". Run run_example.m, ', ...
         'or run RW_G.m and tdRL_G.m before this script.'], outputDir);
end

rwResult = load(rwFile);
tdrlResult = load(tdrlFile);

assert(isequal(rwResult.originalTrials, tdrlResult.originalTrials), ...
    'R-W and tdRL trial coordinates do not match.');
assert(isequal(size(rwResult.observedOriginal), ...
    size(tdrlResult.observedOriginal)), ...
    'R-W and tdRL observed-data dimensions do not match.');
assert(max(abs(rwResult.observedOriginal - ...
    tdrlResult.observedOriginal), [], 'all') < 1e-10, ...
    'R-W and tdRL fits were generated from different observed data.');
assert(isequal(size(rwResult.yfit), size(tdrlResult.yfit)), ...
    'R-W and tdRL prediction dimensions do not match.');

numberOfTrials = 5;
if numel(rwResult.originalTrials) < numberOfTrials
    error('ModelComparison:TooFewTrials', ...
        'At least five measured trials are required for model comparison.');
end

comparisonTrials = rwResult.originalTrials(1:numberOfTrials);
predictionRows = zeros(numberOfTrials, 1);
for trialIndex = 1:numberOfTrials
    [distance, predictionRows(trialIndex)] = min( ...
        abs(rwResult.tspan - comparisonTrials(trialIndex)));
    if distance > 1e-10
        error('ModelComparison:TrialNotFound', ...
            'Measured trial %.6g is absent from the prediction grid.', ...
            comparisonTrials(trialIndex));
    end
end

% Convert freezing percentages to proportions before calculating errors.
observed = rwResult.observedOriginal(1:numberOfTrials, :) / 100;
tdrlPredicted = tdrlResult.yfit(predictionRows, :) / 100;
rwPredicted = rwResult.yfit(predictionRows, :) / 100;
numberOfAnimals = size(observed, 2);

%% SSE, MSE, and RMSE for each animal
tdrlSSE = sum((observed - tdrlPredicted).^2, 1);
rwSSE = sum((observed - rwPredicted).^2, 1);
tdrlMSE = mean((observed - tdrlPredicted).^2, 1);
rwMSE = mean((observed - rwPredicted).^2, 1);
tdrlRMSE = sqrt(tdrlMSE);
rwRMSE = sqrt(rwMSE);

fprintf('\nModel comparison across %d animals and %d measured trials\n', ...
    numberOfAnimals, numberOfTrials);
fprintf('tdRL mean RMSE = %.6f\n', mean(tdrlRMSE));
fprintf('R-W  mean RMSE = %.6f\n', mean(rwRMSE));

%% Normality checks and paired test of animal-level RMSE
testAlpha = 0.05;
[rejectNormalTdrl, normalityPTdrl] = ...
    lillietest(tdrlRMSE, 'Alpha', testAlpha);
[rejectNormalRw, normalityPRw] = ...
    lillietest(rwRMSE, 'Alpha', testAlpha);

degreesOfFreedom = NaN;
confidenceIntervalLower = NaN;
confidenceIntervalUpper = NaN;

if rejectNormalTdrl == 0 && rejectNormalRw == 0
    [~, comparisonP, confidenceInterval, testStats] = ...
        ttest(tdrlRMSE, rwRMSE, 'Alpha', testAlpha);
    testName = "paired t-test";
    testStatistic = testStats.tstat;
    degreesOfFreedom = testStats.df;
    confidenceIntervalLower = confidenceInterval(1);
    confidenceIntervalUpper = confidenceInterval(2);
    fprintf('Both RMSE distributions passed Lilliefors tests.\n');
    fprintf('Paired t-test: t(%d) = %.6f, p = %.6g\n', ...
        degreesOfFreedom, testStatistic, comparisonP);
else
    [comparisonP, ~, testStats] = signrank(tdrlRMSE, rwRMSE, ...
        'Alpha', testAlpha);
    testName = "Wilcoxon signed-rank test";
    testStatistic = testStats.signedrank;
    fprintf('At least one RMSE distribution failed a Lilliefors test.\n');
    fprintf('Wilcoxon signed-rank test: W = %.6f, p = %.6g\n', ...
        testStatistic, comparisonP);
end

%% Save model-comparison outputs
animalId = (1:numberOfAnimals)';
perAnimalTable = table(animalId, tdrlSSE', rwSSE', tdrlMSE', rwMSE', ...
    tdrlRMSE', rwRMSE', 'VariableNames', ...
    {'animal_id', 'tdRL_SSE', 'RW_SSE', 'tdRL_MSE', 'RW_MSE', ...
     'tdRL_RMSE', 'RW_RMSE'});
modelSummaryTable = table(["tdRL"; "R-W"], ...
    [mean(tdrlRMSE); mean(rwRMSE)], ...
    [std(tdrlRMSE); std(rwRMSE)], ...
    repmat(numberOfAnimals, 2, 1), ...
    'VariableNames', {'model', 'mean_RMSE', 'SD_RMSE', 'n_animals'});
normalityTable = table(["tdRL"; "R-W"], ...
    [rejectNormalTdrl; rejectNormalRw], ...
    [normalityPTdrl; normalityPRw], ...
    'VariableNames', {'model', 'reject_normality', 'p_value'});
testTable = table(testName, testStatistic, degreesOfFreedom, comparisonP, ...
    confidenceIntervalLower, confidenceIntervalUpper, testAlpha, ...
    'VariableNames', {'test', 'statistic', 'degrees_of_freedom', ...
     'p_value', 'CI_lower', 'CI_upper', 'alpha'});

writetable(perAnimalTable, ...
    fullfile(outputDir, 'model_comparison_per_animal.csv'));
writetable(modelSummaryTable, ...
    fullfile(outputDir, 'model_comparison_summary.csv'));
writetable(normalityTable, ...
    fullfile(outputDir, 'model_comparison_normality.csv'));
writetable(testTable, ...
    fullfile(outputDir, 'model_comparison_test.csv'));
save(fullfile(outputDir, 'model_comparison_results.mat'), ...
    'comparisonTrials', 'observed', 'tdrlPredicted', 'rwPredicted', ...
    'tdrlSSE', 'rwSSE', 'tdrlMSE', 'rwMSE', 'tdrlRMSE', 'rwRMSE', ...
    'rejectNormalTdrl', 'rejectNormalRw', 'normalityPTdrl', ...
    'normalityPRw', 'testName', 'testStatistic', 'degreesOfFreedom', ...
    'comparisonP', 'testAlpha');

%% Plot mean residuals across measured trials
tdrlResidual = observed - tdrlPredicted;
rwResidual = observed - rwPredicted;
meanTdrlResidual = mean(tdrlResidual, 2);
meanRwResidual = mean(rwResidual, 2);
semTdrlResidual = std(tdrlResidual, 0, 2) / sqrt(numberOfAnimals);
semRwResidual = std(rwResidual, 0, 2) / sqrt(numberOfAnimals);

residualFigure = figure('Color', 'w');
hold on;
errorbar(1:numberOfTrials, meanTdrlResidual, semTdrlResidual, ...
    'r-o', 'LineWidth', 2, 'MarkerSize', 6);
errorbar(1:numberOfTrials, meanRwResidual, semRwResidual, ...
    'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Trial');
ylabel('Mean residual (observed - predicted)');
legend({'tdRL', 'R-W'}, 'Location', 'best');
title('Mean residual across measured trials');
set(gca, 'TickDir', 'out');
exportgraphics(residualFigure, ...
    fullfile(outputDir, 'model_comparison_residuals.png'), ...
    'Resolution', 300);

%% Plot residual distributions
tdrlColor = [146, 71, 131] / 255;
rwColor = [65, 119, 57] / 255;
distributionFigure = figure('Color', 'w');
subplot(1, 2, 1);
histogram(tdrlResidual(:), 'Normalization', 'probability', ...
    'FaceColor', tdrlColor, 'FaceAlpha', 0.7);
xlabel('Residual');
ylabel('Probability');
title('tdRL residual distribution');
subplot(1, 2, 2);
histogram(rwResidual(:), 'Normalization', 'probability', ...
    'FaceColor', rwColor, 'FaceAlpha', 0.7);
xlabel('Residual');
ylabel('Probability');
title('R-W residual distribution');
exportgraphics(distributionFigure, ...
    fullfile(outputDir, 'model_comparison_residual_distributions.png'), ...
    'Resolution', 300);

%% Plot cumulative squared error across measured trials
meanTdrlLoss = mean((observed - tdrlPredicted).^2, 2);
meanRwLoss = mean((observed - rwPredicted).^2, 2);
cumulativeFigure = figure('Color', 'w');
hold on;
plot(1:numberOfTrials, cumsum(meanTdrlLoss), '-o', 'LineWidth', 2);
plot(1:numberOfTrials, cumsum(meanRwLoss), '-s', 'LineWidth', 2);
xlabel('Trial');
ylabel('Cumulative squared error');
legend({'tdRL', 'R-W'}, 'Location', 'northwest');
set(gca, 'TickDir', 'out');
box off;
exportgraphics(cumulativeFigure, ...
    fullfile(outputDir, 'model_comparison_cumulative_error.png'), ...
    'Resolution', 300);
