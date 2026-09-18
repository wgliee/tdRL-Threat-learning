function results = DRL_regression(t, ytrue, mean_resp_eachgroup, varargin)
%% DRL model fitting and response-PDF regression
% This function fits a four-parameter DRL model to behavioral trajectories
%       for individual mice, then tests whether neural responses predict the
%       fitted response probability distribution.
%
% @author:  Mingyang Wei, Shanghai Jiao Tong University, 2026
%
% Inputs:
%       - t: time vector after preprocessing/interpolation (# time points x 1)
%       - ytrue: behavioral freezing trajectories (# time points x # mice)
%       - mean_resp_eachgroup: mean neural responses (# mice x # bins x # groups)
%
% Optional parameters:
%       - 'LowerBound': lower bound for [lambda, k, mu, sigma]
%       - 'UpperBound': upper bound for [lambda, k, mu, sigma]
%       - 'InitialParams': initial guess for [lambda, k, mu, sigma]
%       - 'SampArr': sample indices used to evaluate fitted PDF
%       - 'PooledGroup': response group used for pooled Resp ~ PDF regression
%       - 'LOOGroups': response groups used for leave-one-mouse-out regression
%       - 'ZScoreTrace': whether to z-score neural response per mouse
%       - 'ZScorePDF': whether to z-score fitted PDF per mouse
%       - 'PlotResults': whether to generate figures
%       - 'Verbose': whether to print fitting/regression progress
%
% Outputs:
%       - results: structure containing model fits, regression results,
%                  leave-one-mouse-out prediction results, and statistics
%


%% Optional parameters
p = inputParser;
addParameter(p, 'LowerBound', [0, 0, 0.4, 0]);
addParameter(p, 'UpperBound', [100, 1, 0.4, 2]);
addParameter(p, 'InitialParams', [65, 1, 0.4, 0.68]);
addParameter(p, 'SampArr', 0:4);
addParameter(p, 'PooledGroup', 2);
addParameter(p, 'LOOGroups', [1, 2]);
addParameter(p, 'ZScoreTrace', true);
addParameter(p, 'ZScorePDF', true);
addParameter(p, 'PlotResults', true);
addParameter(p, 'Verbose', true);
parse(p, varargin{:});

lb = p.Results.LowerBound;
ub = p.Results.UpperBound;
p_init = p.Results.InitialParams;
samp_arr = p.Results.SampArr;
pooled_group = p.Results.PooledGroup;
loo_groups = p.Results.LOOGroups;
zscore_trace = p.Results.ZScoreTrace;
zscore_pdf = p.Results.ZScorePDF;
plot_results = p.Results.PlotResults;
verbose = p.Results.Verbose;

%% Input formatting and basic checks
t = t(:);
tspan = t;
numMice = size(ytrue, 2);

if size(ytrue, 1) ~= length(t)
    error('The number of rows in ytrue must match the length of t.');
end

if size(mean_resp_eachgroup, 1) ~= numMice
    error('The first dimension of mean_resp_eachgroup must match the number of mice in ytrue.');
end

%% Fit four DRL parameters for each mouse
% params_all: each row contains [lambda, k, mu, sigma] for one mouse

y0 = ytrue(1,:)';
params_all = zeros(numMice, 4);
sumsq_all = zeros(numMice, 1);
yfit = zeros(length(tspan), numMice);

for mouseIdx = 1:numMice
    if verbose
        fprintf('\n====== Fitting mouse %d ======\n', mouseIdx);
    end

    y_mouse = ytrue(:, mouseIdx);
    y0_mouse = y0(mouseIdx);

    params = optimvar('params', 4, ...
                      'LowerBound', lb, ...
                      'UpperBound', ub);

    % Solve the DRL ODE for this mouse using the current candidate params.
    myfcn_mouse = fcn2optimexpr(@(p_use) RtoODE_matrix(p_use, tspan, y0_mouse, 1), params);

    % Objective: squared error between model trajectory and behavioral trajectory.
    obj_mouse = sum((myfcn_mouse(:) - y_mouse(:)).^2);
    prob_mouse = optimproblem('Objective', obj_mouse);

    r0_mouse.params = p_init;
    [rsol_mouse, ssq_mouse] = solve(prob_mouse, r0_mouse);

    params_all(mouseIdx, :) = rsol_mouse.params(:)';
    sumsq_all(mouseIdx) = ssq_mouse;

    if verbose
        fprintf('Mouse %d: params = [%.4f  %.4f  %.4f  %.4f], SSE = %.4f\n', ...
            mouseIdx, params_all(mouseIdx, :), ssq_mouse);
    end

    p_fit_mouse = rsol_mouse.params;
    sol_mouse = ode45(@(t_use, y_use) drldiffun(t_use, y_use, p_fit_mouse), tspan, y0_mouse);
    yfit(:, mouseIdx) = deval(sol_mouse, tspan);
end

%% Prepare fitted PDF for each mouse
pdf_eachmice = zeros(numMice, length(samp_arr));
for micenum = 1:numMice
    param_temp = params_all(micenum, :);
    pdf_eachmice(micenum, :) = myNormalPDF(samp_arr, param_temp(3), param_temp(4));
end

%% Pooled response-PDF regression
% Resp = beta0 + beta1 * PDF

resp_use = mean_resp_eachgroup(:, :, pooled_group);
resp_use = zscore(resp_use, 0, 2);

PDF_all = reshape(pdf_eachmice', [], 1);
Resp_all = reshape(resp_use', [], 1);

tbl = table(PDF_all, Resp_all, 'VariableNames', {'PDF', 'Resp'});
mdl = fitlm(tbl, 'Resp ~ PDF');

if verbose
    disp(mdl);
end

beta0 = mdl.Coefficients{'(Intercept)', 'Estimate'};
beta1 = mdl.Coefficients{'PDF', 'Estimate'};
p_slope = mdl.Coefficients{'PDF', 'pValue'};
R2 = mdl.Rsquared.Ordinary;

if verbose
    fprintf('Resp = %.3f + %.3f * PDF,  slope p=%.4g,  R^2=%.3f\n', ...
        beta0, beta1, p_slope, R2);
end

if plot_results
    figure;
    scatter(PDF_all, Resp_all, 'filled'); hold on;

    xfit = linspace(min(PDF_all), max(PDF_all), 100)';
    yfit_all = predict(mdl, table(xfit, 'VariableNames', {'PDF'}));
    plot(xfit, yfit_all, 'LineWidth', 2);

    xlabel('PDF value');
    ylabel('Neural response');
    title('All mice pooled: linear regression Resp ~ PDF');
    grid on;

    txt = sprintf('slope p = %.3g\nR^2 = %.3f', p_slope, R2);
    xpos = min(PDF_all) + 0.05 * (max(PDF_all) - min(PDF_all));
    ypos = max(Resp_all) - 0.1 * (max(Resp_all) - min(Resp_all));
    text(xpos, ypos, txt, 'FontSize', 10, 'BackgroundColor', 'w');
    hold off;
end

%% Response-PDF leave-one-mouse-out regression
% PDF = beta0 + beta1 * Resp

pdf_mu_eachmice = mean(pdf_eachmice, 2);
pdf_std_eachmice = std(pdf_eachmice, 0, 2);

if zscore_pdf
    pdf_eachmice_use = (pdf_eachmice - pdf_mu_eachmice) ./ pdf_std_eachmice;
else
    pdf_eachmice_use = pdf_eachmice;
end

% Average 5 bins into early and late bins.
pdf_earlylate_eachmice = [mean(pdf_eachmice_use(:, 1:2), 2), ...
                          mean(pdf_eachmice_use(:, 3:5), 2)];
samp_arr_earlylate = 0:1;

loo_results = struct();
for ii = 1:length(loo_groups)
    group_id = loo_groups(ii);
    group_name = sprintf('g%d', group_id);

    resp_use = mean_resp_eachgroup(:, :, group_id);
    if zscore_trace
        resp_use = zscore(resp_use, 0, 2);
    end

    resp_use_earlylate = [mean(resp_use(:, 1:2), 2), ...
                          mean(resp_use(:, 3:5), 2)];

    loo_results.(group_name) = run_loo_regression(resp_use, resp_use_earlylate, ...
        pdf_earlylate_eachmice, samp_arr_earlylate, group_name, verbose, plot_results);
end

%% Compare LOO error with dummy-regressor baseline
% Same final statistics as in the original script:
%       [~,p] = ttest(mae_median_bl_all_40(:), abs_err_g1(:))
%       [~,p] = ttest(abs_err_g2(:), mae_median_bl_all_40(:))
%
% The p values are stored with explicit names in results.stats.

stats = struct();
if isfield(loo_results, 'g1')
    mae_median_bl_all_40 = loo_results.g1.mae_median_bl_all;
    abs_err_g1 = loo_results.g1.abs_err;

    [~, p_g1_vs_median_baseline] = ttest(mae_median_bl_all_40(:), abs_err_g1(:));

    stats.mae_median_bl_all_40 = mae_median_bl_all_40;
    stats.abs_err_g1 = abs_err_g1;
    stats.p_g1_vs_median_baseline = p_g1_vs_median_baseline;

    if verbose
        fprintf('g1 abs error vs median baseline: p = %.4g\n', p_g1_vs_median_baseline);
    end

    if isfield(loo_results, 'g2')
        abs_err_g2 = loo_results.g2.abs_err;
        [~, p_g2_vs_median_baseline] = ttest(abs_err_g2(:), mae_median_bl_all_40(:));

        stats.abs_err_g2 = abs_err_g2;
        stats.p_g2_vs_median_baseline = p_g2_vs_median_baseline;

        if verbose
            fprintf('g2 abs error vs median baseline: p = %.4g\n', p_g2_vs_median_baseline);
        end
    end
end

%% Collect outputs
results = struct();
results.t = t;
results.ytrue = ytrue;
results.y0 = y0;
results.params_all = params_all;
results.sumsq_all = sumsq_all;
results.yfit = yfit;
results.pdf_eachmice = pdf_eachmice;
results.pdf_earlylate_eachmice = pdf_earlylate_eachmice;

results.pooled = struct();
results.pooled.group = pooled_group;
results.pooled.mdl = mdl;
results.pooled.beta0 = beta0;
results.pooled.beta1 = beta1;
results.pooled.p_slope = p_slope;
results.pooled.R2 = R2;
results.pooled.PDF_all = PDF_all;
results.pooled.Resp_all = Resp_all;

results.loo = loo_results;
results.stats = stats;

% Compatibility fields matching the original workspace variable names.
if isfield(stats, 'mae_median_bl_all_40')
    results.mae_median_bl_all_40 = stats.mae_median_bl_all_40;
end
if isfield(stats, 'abs_err_g1')
    results.abs_err_g1 = stats.abs_err_g1;
end
if isfield(stats, 'abs_err_g2')
    results.abs_err_g2 = stats.abs_err_g2;
end

results.settings = p.Results;
end

function out = run_loo_regression(resp_use, resp_use_earlylate, pdf_earlylate_eachmice, ...
    samp_arr_earlylate, group_name, verbose, plot_results)
%% Leave-one-mouse-out linear regression: PDF = beta0 + beta1 * Resp

[nMice, nBins] = size(pdf_earlylate_eachmice);

beta_cv = zeros(nMice, 2);
R_test_each = zeros(nMice, 1);
R2_test_each = zeros(nMice, 1);
abs_err = zeros(nMice, nBins);
sse_each = zeros(nMice, 1);
pdf_pred_eachmice = zeros(nMice, nBins);
pdf_pred_Full = zeros(nMice, size(resp_use, 2));
mae_median_bl = zeros(nMice, 1);
mae_mean_bl = zeros(nMice, 1);
mae_median_bl_all = zeros(nMice, nBins);

for m = 1:nMice
    if verbose
        fprintf('\n=== LOO %s: holding out mouse %d as test ===\n', group_name, m);
    end

    train_idx = setdiff(1:nMice, m);

    Xtrain = resp_use_earlylate(train_idx, :);
    Ytrain = pdf_earlylate_eachmice(train_idx, :);

    Xtrain = Xtrain(:);
    Ytrain = Ytrain(:);

    tblTrain = table(Xtrain, Ytrain, 'VariableNames', {'Resp', 'PDF'});
    mdl = fitlm(tblTrain, 'PDF ~ Resp');

    if verbose
        disp(mdl);
    end

    beta0 = mdl.Coefficients{'(Intercept)', 'Estimate'};
    beta1 = mdl.Coefficients{'Resp', 'Estimate'};
    beta_cv(m, :) = [beta0, beta1];

    Xtest = resp_use_earlylate(m, :)';
    Ytest = pdf_earlylate_eachmice(m, :)';

    tblTest = table(Xtest, 'VariableNames', {'Resp'});
    Ypred = predict(mdl, tblTest);
    pdf_pred_eachmice(m, :) = Ypred';

    Resp_test_full = resp_use(m, :)';
    tblTestFull = table(Resp_test_full, 'VariableNames', {'Resp'});
    Resp_pred = predict(mdl, tblTestFull);
    pdf_pred_Full(m, :) = Resp_pred';

    R_test_each(m) = corr(Ypred, Ytest);
    abs_err(m, :) = sqrt((Ytest - Ypred).^2)';

    SSE = sum((Ytest - Ypred).^2);
    SST = sum((Ytest - mean(Ytest)).^2);
    sse_each(m) = SSE;
    R2_test_each(m) = 1 - SSE / SST;

    mae_median_bl(m) = mean(abs(Ytest - median(Ytrain)), 1);
    mae_median_bl_all(m, :) = abs(Ytest - median(Ytrain))';
    mae_mean_bl(m) = mean(abs(Ytest - mean(Ytrain)), 1);

    if verbose
        fprintf('Mouse %d: test R = %.3f, R^2 = %.3f\n', ...
            m, R_test_each(m), R2_test_each(m));
    end
end

if verbose
    fprintf('\n=== LOO %s summary ===\n', group_name);
    fprintf('Mean test R   = %.3f\n', mean(R_test_each));
    fprintf('Mean test R^2 = %.3f\n', mean(R2_test_each));
end

if plot_results
    figure;
    tiledlayout(2, ceil(nMice/2), 'TileSpacing', 'compact', 'Padding', 'compact');

    for m = 1:nMice
        nexttile;
        plot(samp_arr_earlylate, pdf_earlylate_eachmice(m, :), 'o-', 'DisplayName', 'True'); hold on;
        plot(samp_arr_earlylate, pdf_pred_eachmice(m, :), '-.', 'DisplayName', 'Predicted');
        title(sprintf('Mouse %d  (R=%.2f, R^2=%.2f)', ...
            m, R_test_each(m), R2_test_each(m)));
        xlabel('Sample index');
        ylabel(['PDF ', group_name]);
        grid on;
        if m == 1
            legend('Location', 'best');
        end
    end
end

out = struct();
out.beta_cv = beta_cv;
out.R_test_each = R_test_each;
out.R2_test_each = R2_test_each;
out.abs_err = abs_err;
out.sse_each = sse_each;
out.mse = sse_each;
out.pdf_pred_eachmice = pdf_pred_eachmice;
out.pdf_pred_Full = pdf_pred_Full;
out.mae_median_bl = mae_median_bl;
out.mae_mean_bl = mae_mean_bl;
out.mae_median_bl_all = mae_median_bl_all;
out.mean_test_R = mean(R_test_each);
out.mean_test_R2 = mean(R2_test_each);
end
