%% Fit the Rescorla-Wagner (R-W) model
% This script can be run directly. For the complete workflow, run
% run_example.m or call run_analysis(inputFile, outputDir).

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'function'));

if ~exist('inputFile', 'var') || isempty(inputFile)
    inputFile = fullfile(scriptDir, 'Example data.xlsx');
end
if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = fullfile(scriptDir, 'results');
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

interpolationStep = 0.2;
[originalTrials, observedOriginal, tspan, ytrue] = ...
    prepare_behavior_data(inputFile, interpolationStep);

numMice = size(ytrue, 2);
y0 = ytrue(1, :)';
yMean = mean(ytrue, 2);
ySE = std(ytrue, 0, 2) / sqrt(numMice);

%% Define the optimization problem
% params(1) = lambda: asymptotic freezing level, range [0, 100]
% params(2) = alpha: learning rate, range [0, 1]
params = optimvar('params', 2, ...
    'LowerBound', [0, 0], ...
    'UpperBound', [100, 1]);

initialParameters = [50, 0.1];
modelExpression = fcn2optimexpr( ...
    @(p) rwRtoODE_matrix(p, tspan, y0, numMice), params);
objective = sum((mean(modelExpression, 2) - yMean).^2);
problem = optimproblem('Objective', objective);
initialPoint.params = initialParameters;

%% Fit and simulate each animal using the shared group-level parameters
[solution, sumSquaredError, exitFlag, solverOutput] = ...
    solve(problem, initialPoint);
pFit = solution.params;
exitFlagValue = double(exitFlag);
yfit = rwRtoODE_matrix(pFit, tspan, y0, numMice);

fprintf('\nR-W model parameters\n');
fprintf('lambda = %.6f\n', pFit(1));
fprintf('alpha  = %.6f\n', pFit(2));
fprintf('group-level SSE = %.6f\n\n', sumSquaredError);

%% Save machine-readable outputs
predictionTable = make_behavior_table(tspan, yfit);
observedTable = make_behavior_table(tspan, ytrue);
parameterTable = table(pFit(1), pFit(2), sumSquaredError, exitFlagValue, ...
    'VariableNames', {'lambda', 'alpha', 'group_level_SSE', 'exit_flag'});

writetable(predictionTable, fullfile(outputDir, 'RW_predictions.csv'));
writetable(observedTable, fullfile(outputDir, 'observed_interpolated.csv'));
writetable(parameterTable, fullfile(outputDir, 'RW_fit_parameters.csv'));
save(fullfile(outputDir, 'RW_fit_results.mat'), ...
    'inputFile', 'originalTrials', 'observedOriginal', 'tspan', 'ytrue', ...
    'yfit', 'yMean', 'ySE', 'pFit', 'sumSquaredError', 'exitFlagValue', ...
    'solverOutput', 'interpolationStep');

%% Plot observed and fitted group means
meanCurveFit = mean(yfit, 2);
fitFigure = figure('Color', 'w');
fill([tspan; flipud(tspan)], ...
    [yMean + ySE; flipud(yMean - ySE)], ...
    [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
observedLine = plot(tspan, yMean, 'LineWidth', 1.5, ...
    'Color', [0.55, 0.55, 0.55]);
fittedLine = plot(tspan, meanCurveFit, 'LineWidth', 2, ...
    'Color', [0.15, 0.15, 0.15]);

xlabel('Trial coordinate');
ylabel('Freezing (%)');
ylim([0, 100]);
xlim([tspan(1), tspan(end)]);
xticks(originalTrials);
title('R-W model: observed and fitted group means');
legend([observedLine, fittedLine], {'Observed mean', 'Fitted mean'}, ...
    'Location', 'best');
set(gca, 'TickDir', 'out');
exportgraphics(fitFigure, fullfile(outputDir, 'RW_fit.png'), ...
    'Resolution', 300);
