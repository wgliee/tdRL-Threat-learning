%% Fit the trial-wise discount reinforcement learning (tdRL) model
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
% params(3) = mu: discount-function center on model coordinates [0, 6]
% params(4) = sigma: discount-function width, range (0, 20]
minimumSigma = 1e-6;
params = optimvar('params', 4, ...
    'LowerBound', [0, 0, 0, minimumSigma], ...
    'UpperBound', [100, 1, 6, 20]);

initialParameters = [50, 0.1, 0.1, 0.4];
modelExpression = fcn2optimexpr( ...
    @(p) RtoODE_matrix(p, tspan, y0, numMice), params);
objective = sum((mean(modelExpression, 2) - yMean).^2);
problem = optimproblem('Objective', objective);
initialPoint.params = initialParameters;

%% Fit and simulate each animal using the shared group-level parameters
[solution, sumSquaredError, exitFlag, solverOutput] = ...
    solve(problem, initialPoint);
pFit = solution.params;
exitFlagValue = double(exitFlag);
yfit = RtoODE_matrix(pFit, tspan, y0, numMice);

% The spreadsheet uses coordinates 0-4 for displayed trials 1-5.
muTrial = pFit(3) + 1;
fprintf('\ntdRL model parameters\n');
fprintf('lambda   = %.6f\n', pFit(1));
fprintf('alpha    = %.6f\n', pFit(2));
fprintf('mu       = %.6f (model coordinate)\n', pFit(3));
fprintf('mu_trial = %.6f (displayed trial number)\n', muTrial);
fprintf('sigma    = %.6f\n', pFit(4));
fprintf('group-level SSE = %.6f\n\n', sumSquaredError);

%% Save machine-readable outputs
predictionTable = make_behavior_table(tspan, yfit);
observedTable = make_behavior_table(tspan, ytrue);
parameterTable = table(pFit(1), pFit(2), pFit(3), muTrial, pFit(4), ...
    sumSquaredError, exitFlagValue, 'VariableNames', ...
    {'lambda', 'alpha', 'mu_model_coordinate', 'mu_trial', 'sigma', ...
     'group_level_SSE', 'exit_flag'});

writetable(predictionTable, fullfile(outputDir, 'tdRL_predictions.csv'));
writetable(observedTable, fullfile(outputDir, 'observed_interpolated.csv'));
writetable(parameterTable, fullfile(outputDir, 'tdRL_fit_parameters.csv'));
save(fullfile(outputDir, 'tdRL_fit_results.mat'), ...
    'inputFile', 'originalTrials', 'observedOriginal', 'tspan', 'ytrue', ...
    'yfit', 'yMean', 'ySE', 'pFit', 'muTrial', 'sumSquaredError', ...
    'exitFlagValue', 'solverOutput', 'interpolationStep');

%% Plot observed and fitted group means
meanCurveFit = mean(yfit, 2);
fitFigure = figure('Color', 'w');
fill([tspan; flipud(tspan)], ...
    [yMean + ySE; flipud(yMean - ySE)], ...
    [0.6, 0.8, 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
observedLine = plot(tspan, yMean, 'LineWidth', 1.5, ...
    'Color', [0.35, 0.55, 0.75]);
fittedLine = plot(tspan, meanCurveFit, 'LineWidth', 2, ...
    'Color', [0, 0.25, 0.75]);

xlabel('Trial coordinate');
ylabel('Freezing (%)');
ylim([0, 100]);
xlim([tspan(1), tspan(end)]);
xticks(originalTrials);
title('tdRL model: observed and fitted group means');
legend([observedLine, fittedLine], {'Observed mean', 'Fitted mean'}, ...
    'Location', 'best');
set(gca, 'TickDir', 'out');
exportgraphics(fitFigure, fullfile(outputDir, 'tdRL_fit.png'), ...
    'Resolution', 300);
