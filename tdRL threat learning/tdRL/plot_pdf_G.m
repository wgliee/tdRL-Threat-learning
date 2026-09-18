%% Plot the normalized Gaussian discount function used by the tdRL model
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'function'));

displayedTrial = 1:0.01:5;
muTrial = 1.40;
sigma = 0.6759;

normalDensity = myNormalPDF(displayedTrial, muTrial, sigma);
discountWeight = normalDensity / max(normalDensity);

discountFigure = figure('Color', 'w');
plot(displayedTrial, discountWeight, 'LineWidth', 1.5, ...
    'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Trial');
ylabel('Discount weight, \phi(t)');
title('Trial-wise discount function');
ylim([0, 1.05]);
grid on;
set(gca, 'TickDir', 'out');

outputDir = fullfile(scriptDir, 'results');
if ~isfolder(outputDir)
    mkdir(outputDir);
end
exportgraphics(discountFigure, ...
    fullfile(outputDir, 'tdRL_discount_function.png'), ...
    'Resolution', 300);
