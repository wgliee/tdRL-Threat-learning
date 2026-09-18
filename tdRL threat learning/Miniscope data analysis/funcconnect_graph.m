function results = funcconnect_graph(traceeach_zs, freezing_each, clustidx_each, varargin)
%% Functional connectivity and graph analysis during CS trials
% This function computes normalized functional connectivity matrices and
% graph-based degree centrality from z-scored neural activity traces.
%
% @author:  Mingyang Wei, Shanghai Jiao Tong University, 2026
%
% Inputs:
%       traceeach_zs   - cell array, traceeach_zs{mouse}: neurons x frames
%                        z-scored neural activity traces for each mouse
%       freezing_each  - cell array, freezing_each{mouse}: behavioral
%                        freezing vector for each mouse
%       clustidx_each  - cell array, clustidx_each{mouse}: cluster labels
%                        for neurons from each mouse
%
% Optional name-value inputs:
%       'DsTimes'           - downsampling factor. Default: 1
%       'SdThreshold'       - threshold for functional connection after
%                             normalization to shuffled distribution.
%                             Default: 1.4
%       'GroupInPaper'      - cell groups included in the analysis.
%                             Default: [1 4]
%       'ToneIncludeLength' - number of frames included after CS onset.
%                             Default: 100
%       'UseAbsAdj'         - if true, threshold abs(normalized FC).
%                             Default: false
%       'PlotResults'       - if true, plot summary figures. Default: true
%       'Verbose'           - if true, print progress. Default: true
%
% Outputs:
%       results - structure containing normalized connectivity matrices,
%                 adjacency matrices, graph metrics, degree centrality,
%                 and summary statistics.
%


%% Input parser
p = inputParser;
addRequired(p, 'traceeach_zs', @(x) iscell(x));
addRequired(p, 'freezing_each', @(x) iscell(x));
addRequired(p, 'clustidx_each', @(x) iscell(x));
addParameter(p, 'DsTimes', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'SdThreshold', 1.4, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'GroupInPaper', [1 4], @(x) isnumeric(x) || isempty(x));
addParameter(p, 'ToneIncludeLength', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'UseAbsAdj', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'PlotResults', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, traceeach_zs, freezing_each, clustidx_each, varargin{:});

DsTimes = p.Results.DsTimes;
sd_threshold = p.Results.SdThreshold;
group_in_paper = p.Results.GroupInPaper;
tone_include_length = p.Results.ToneIncludeLength;
use_abs_adj = logical(p.Results.UseAbsAdj);
plot_results = logical(p.Results.PlotResults);
verbose = logical(p.Results.Verbose);

n = numel(traceeach_zs);

%% CS trial indices
% CS window fixed from the original script.
tone_start = ([180 330 450 630 720] * 10) ./ DsTimes;
tone_start = [tone_start, ([180 260 350 410 490 570 630 680 730 780 870 930] * 10 + 7600) ./ DsTimes];
tone_start = round(tone_start);
numTones = numel(tone_start);

% Kept from the original code, although shock data are not used for the
% functional connectivity calculation below.
shock_start = ([180 330 450 630 720] * 10 + 300) ./ DsTimes;
shock_start = round(shock_start);
shock_include_length = round(40 ./ DsTimes);

%% Functional connectivity for CS trials
trail_freezing = [];
adjacency_mat_allmice = cell(n, 1);
adjacency_mat_sorted_all = cell(n, 1);
norm_connectivity_mat_allmice = cell(n, 1);
norm_sorted_connectivity_mat_all = cell(n, 1);
labelsXY_eachmice = cell(n, 1);
endgroup_idx_all = [];
connect = [];
noise_connect = [];
corr_temp = [];
mean_freeze = [];

for micenum = 1:n
    if verbose
        fprintf('\nNow processing micenum: %d', micenum);
    end

    if DsTimes == 1
        micetrace = traceeach_zs{micenum};
        micefreeze = freezing_each{micenum};
    else
        micetrace = downsample_trace(traceeach_zs{micenum}, DsTimes);
        micefreeze = downsample_trace(freezing_each{micenum}', DsTimes)';
    end

    micefreeze = micefreeze(:)';
    curmiceidx = clustidx_each{micenum};
    curmiceidx = curmiceidx(:);

    % Select neuronal subgroups used in the original analysis.
    if ~isempty(group_in_paper)
        keep_cell = ismember(curmiceidx, group_in_paper);
        micetrace = micetrace(keep_cell, :);
        curmiceidx = curmiceidx(keep_cell);
    end

    tone_dataset = [];
    pretone_dataset = [];
    shock_dataset = [];
    freeze_dataset = [];

    % Prepare CS dataset: neurons x frames x trials.
    for tones = 1:numTones
        tone_dataset = cat(3, tone_dataset, micetrace(:, tone_start(tones):tone_start(tones)+tone_include_length-1));
        pretone_dataset = cat(3, pretone_dataset, micetrace(:, tone_start(tones)-tone_include_length:tone_start(tones)-1));
        freeze_dataset = [freeze_dataset, micefreeze(tone_start(tones):tone_start(tones)+tone_include_length-1)];
    end

    mean_freeze(micenum, :) = mean(freeze_dataset, 1);
    trail_freezing(micenum, :) = mean(freeze_dataset, 1);

    for shocks = 1:numel(shock_start)
        shock_dataset = cat(3, shock_dataset, micetrace(:, shock_start(shocks):shock_start(shocks)+shock_include_length-1));
    end

    % Correlation of each trial response with every other trial response.
    for i = 1:size(tone_dataset, 3)
        for j = 1:size(tone_dataset, 3)
            corr_temp(i, j, micenum) = mean(diag(corr(tone_dataset(:, :, i), tone_dataset(:, :, j))));
        end
    end

    % Functional connectivity.
    mean_tone = mean(tone_dataset, 3)';
    stacked_mean_tone = repmat(mean(tone_dataset, 3), [1 1 size(tone_dataset, 3)]);

    % Sort neurons by selected functional groups for visualization and
    % in-/out-group connectivity analysis.
    labelsXY = repmat({""}, size(tone_dataset, 1), 1);
    sort_idx = [];
    endgroup_idx = [];
    for i = 1:length(group_in_paper)
        sort_idx = [sort_idx; find(curmiceidx == group_in_paper(i))];
        endgroup_idx(i, 1) = length(sort_idx);
        labelsXY{endgroup_idx(i, 1)} = string(endgroup_idx(i, 1));
    end
    labelsXY_eachmice{micenum} = labelsXY;
    endgroup_idx_all = [endgroup_idx_all, [0; endgroup_idx]];

    if isempty(sort_idx)
        sorted_tone_dataset = tone_dataset;
    else
        sorted_tone_dataset = tone_dataset(sort_idx, :, :);
    end

    % Compute mean and SD distribution for normalization.
    [meanMatrix, stdMatrix] = generateMeanStdMatrices(tone_dataset);
    [meandenoisedMatrix, stddenoisedMatrix] = generateMeanStdMatrices(tone_dataset - stacked_mean_tone);
    [meansortedMatrix, stdsortedMatrix] = generateMeanStdMatrices(sorted_tone_dataset);

    adjacency_mat_all = [];
    norm_connectivity_mat_all = [];

    for tones = 1:size(tone_dataset, 3)
        nodes = size(meanMatrix, 1);
        elements_mask = ~logical(eye(nodes));

        % Normalization factors from shuffled distributions.
        meantemp = meanMatrix(:, :, tones);
        stdtemp = stdMatrix(:, :, tones);
        meandenoisedtemp = meandenoisedMatrix(:, :, tones);
        stddenoisedtemp = stddenoisedMatrix(:, :, tones);
        meansorttemp = meansortedMatrix(:, :, tones);
        stdsorttemp = stdsortedMatrix(:, :, tones);

        temp = tone_dataset(:, :, tones)';
        noise_temp = temp - mean_tone;
        sorted_temp = sorted_tone_dataset(:, :, tones)';

        connectivity_mat = corr(temp, temp);
        noise_mat = corr(noise_temp, noise_temp);
        sorted_connectivity_mat = corr(sorted_temp, sorted_temp);

        % Normalize actual correlation relative to shuffled distribution.
        norm_connectivity_mat = (connectivity_mat - meantemp) ./ stdtemp;
        norm_noise_mat = (noise_mat - meandenoisedtemp) ./ stddenoisedtemp;
        norm_sorted_connectivity_mat = (sorted_connectivity_mat - meansorttemp) ./ stdsorttemp;

        norm_connectivity_mat(~elements_mask) = 1;
        norm_noise_mat(~elements_mask) = 1;
        norm_sorted_connectivity_mat(logical(eye(size(norm_sorted_connectivity_mat, 1)))) = 1;

        norm_sorted_connectivity_mat_all{micenum}(:, :, tones) = norm_sorted_connectivity_mat;

        conn_distribution = norm_connectivity_mat(elements_mask);
        noise_distribution = norm_noise_mat(elements_mask);
        connect(micenum, tones) = mean(conn_distribution);
        noise_connect(micenum, tones) = mean(noise_distribution);

        % Sorted adjacency matrix.
        acyclic_cn_mat = norm_sorted_connectivity_mat;
        acyclic_cn_mat(logical(eye(size(acyclic_cn_mat, 1)))) = 0;
        adjacency_mat = acyclic_cn_mat;

        if ~use_abs_adj
            connected_idx = adjacency_mat > sd_threshold;
        else
            connected_idx = abs(adjacency_mat) > sd_threshold;
        end

        adjacency_mat(connected_idx) = 1;
        adjacency_mat(~connected_idx) = 0;
        adjacency_mat_sorted_all{micenum}(:, :, tones) = adjacency_mat;

        % Unsorted adjacency matrix used for graph analysis.
        acyclic_cn_mat = norm_connectivity_mat;
        acyclic_cn_mat(logical(eye(nodes))) = 0;
        adjacency_mat = acyclic_cn_mat;

        if ~use_abs_adj
            connected_idx = adjacency_mat > sd_threshold;
        else
            connected_idx = abs(adjacency_mat) > sd_threshold;
        end

        adjacency_mat(connected_idx) = 1;
        adjacency_mat(~connected_idx) = 0;
        adjacency_mat_all(:, :, tones) = adjacency_mat;
        norm_connectivity_mat_all(:, :, tones) = norm_connectivity_mat;
    end

    adjacency_mat_allmice{micenum} = adjacency_mat_all;
    norm_connectivity_mat_allmice{micenum} = norm_connectivity_mat_all;
end

%% In-group and out-group connectivity
meantriu_tones = [];
mean_outgroup_tones = [];

for micenum = 1:n
    for tones = 1:numTones
        mat_use = norm_sorted_connectivity_mat_all{micenum}(:, :, tones);
        endgroup_idx = endgroup_idx_all(:, micenum);

        for groupcount = 1:length(group_in_paper)
            groups = group_in_paper(groupcount);
            current_indices = endgroup_idx(groupcount)+1:endgroup_idx(groupcount+1);

            submat = mat_use(current_indices, current_indices);
            uppertri = triu(submat, 1);

            if nnz(uppertri) == 0
                meantriu_tones(micenum, tones, groups) = 0;
            else
                meantriu_tones(micenum, tones, groups) = sum(uppertri, "all") / nnz(uppertri);
            end

            outgroup_mat = mat_use(current_indices, :);
            outgroup_mat(:, current_indices) = 0;
            total_outgroup = sum(outgroup_mat, 'all');
            count_outgroup = nnz(outgroup_mat);

            if count_outgroup == 0
                mean_outgroup_tones(micenum, tones, groups) = 0;
            else
                mean_outgroup_tones(micenum, tones, groups) = total_outgroup / count_outgroup;
            end
        end
    end
end

g1_within = mean(meantriu_tones(:, :, 1), 3);
g4_within = mean(meantriu_tones(:, :, 4), 3);
g1_zs = zscore(g1_within')';
g4_zs = zscore(g4_within')';

%% Graph analysis
centrallity_all_mice = cell(n, 1);
degree_all_together = [];
degree_all_mice = cell(n, 1);
centrality_each_group = NaN(n, 4);
numComponents = NaN(n, numTones);
charPathLength = NaN(n, numTones);
networkDiameter = NaN(n, numTones);
globalEfficiency = NaN(n, numTones);
density = NaN(n, numTones);
meanClusteringCoeff = NaN(n, numTones);

for micenum = 1:n
    centrality_all = [];
    degree_all = [];

    curmiceidx = clustidx_each{micenum};
    curmiceidx = curmiceidx(:);
    if ~isempty(group_in_paper)
        curmiceidx(~ismember(curmiceidx, group_in_paper), :) = [];
    end

    cur_cell_num = length(curmiceidx);
    adjacency_mat_all = adjacency_mat_allmice{micenum};

    for tones = 1:numTones
        adjacency_temp = adjacency_mat_all(:, :, tones);
        G = graph(adjacency_temp);

        [bin, ~] = conncomp(G);
        numComponents(micenum, tones) = max(bin) ./ length(adjacency_temp);

        D = distances(G);
        D_upper = triu(D, 1);
        charPathLength(micenum, tones) = mean(D_upper(D_upper ~= 0 & ~isinf(D_upper)));
        networkDiameter(micenum, tones) = max(D(~isinf(D)));

        finiteDist = D(~isinf(D) & D ~= 0);
        globalEfficiency(micenum, tones) = mean(1 ./ finiteDist);
        density(micenum, tones) = numedges(G) / (numnodes(G) * (numnodes(G) - 1) / 2);

        centrality_temp = centrality(G, 'degree')';
        centrality_all(:, tones) = centrality_temp;
        degree_all(:, tones) = centrality_temp;

        clustering_coeffs = calculate_clustering_coeff(adjacency_temp);
        meanClusteringCoeff(micenum, tones) = mean(clustering_coeffs);
    end

    centrallity_all_mice{micenum} = centrality_all;
    degree_all_together = [degree_all_together; degree_all .* 100 ./ cur_cell_num];
    degree_all_mice{micenum} = degree_all .* 100 ./ cur_cell_num;

    norm_centrality_all = zscore(centrality_all);
    for g = 1:4
        centrality_each_group(micenum, g) = mean(norm_centrality_all(curmiceidx == g, :), [1, 2], 'omitnan');
    end
end

%% Degree analysis by group
centrality_each_group_rowzs = NaN(n, numTones, 4);

for micenum = 1:n
    curmiceidx = clustidx_each{micenum};
    curmiceidx = curmiceidx(:);
    if ~isempty(group_in_paper)
        curmiceidx(~ismember(curmiceidx, group_in_paper), :) = [];
    end

    centrality_all = centrallity_all_mice{micenum};
    centrality_all = zscore(centrality_all);

    for tones = 1:numTones
        for g = 1:4
            centrality_each_group_rowzs(micenum, tones, g) = mean(centrality_all(curmiceidx == g, tones), 1, 'omitnan');
        end
    end
end

g1_centrality = mean(centrality_each_group_rowzs(:, :, 1), 3, 'omitnan');
g4_centrality = mean(centrality_each_group_rowzs(:, :, 4), 3, 'omitnan');

[~, p_degree_early] = ttest(mean(g1_centrality(:, 1:2), 2, 'omitnan'), mean(g4_centrality(:, 1:2), 2, 'omitnan'));
[~, p_degree_late] = ttest(mean(g1_centrality(:, 3:5), 2, 'omitnan'), mean(g4_centrality(:, 3:5), 2, 'omitnan'));

%% Optional plotting
if plot_results
    figure;
    rowdata_columnumber_plot(connect);
    title("connectivity")

    figure;
    rowdata_columnumber_plot(noise_connect);
    title("denoised connectivity")

    figure;
    h = heatmap(mean(corr_temp, 3));
    h.GridVisible = 0;
    h.Colormap = jet;
    h.ColorLimits = [0, 0.1];

    figure;
    subplot(211)
    hold on;
    rowdata_columnumber_plot(g1_zs);
    rowdata_columnumber_plot(g4_zs);
    hold off;
    legend("1", "4");
    title('In-group Connectivity');

    usetest = mean_outgroup_tones(:, :, 1);
    usetest = zscore(usetest')';
    subplot(212);
    hold on;
    rowdata_columnumber_plot(mean(usetest(:, :, 1), 3));
    hold off;
    legend("1");
    title('Out-group Connectivity');

    figure;
    hold on
    rowdata_columnumber_plot(g1_centrality);
    rowdata_columnumber_plot(g4_centrality);
    legend("g1", "g4")
    hold off
end

%% Collect outputs
results = struct();
results.params.DsTimes = DsTimes;
results.params.SdThreshold = sd_threshold;
results.params.GroupInPaper = group_in_paper;
results.params.ToneIncludeLength = tone_include_length;
results.params.UseAbsAdj = use_abs_adj;
results.params.ToneStart = tone_start;
results.params.Window = "CS";

results.trail_freezing = trail_freezing;
results.mean_freeze = mean_freeze;
results.corr_temp = corr_temp;

results.connect = connect;
results.noise_connect = noise_connect;
results.norm_connectivity_mat_allmice = norm_connectivity_mat_allmice;
results.norm_sorted_connectivity_mat_all = norm_sorted_connectivity_mat_all;
results.adjacency_mat_allmice = adjacency_mat_allmice;
results.adjacency_mat_sorted_all = adjacency_mat_sorted_all;
results.labelsXY_eachmice = labelsXY_eachmice;
results.endgroup_idx_all = endgroup_idx_all;

results.group_connectivity.meantriu_tones = meantriu_tones;
results.group_connectivity.mean_outgroup_tones = mean_outgroup_tones;
results.group_connectivity.g1_within = g1_within;
results.group_connectivity.g4_within = g4_within;
results.group_connectivity.g1_zs = g1_zs;
results.group_connectivity.g4_zs = g4_zs;

results.graph.centrallity_all_mice = centrallity_all_mice;
results.graph.degree_all_together = degree_all_together;
results.graph.degree_all_mice = degree_all_mice;
results.graph.centrality_each_group = centrality_each_group;
results.graph.centrality_each_group_rowzs = centrality_each_group_rowzs;
results.graph.g1_centrality = g1_centrality;
results.graph.g4_centrality = g4_centrality;
results.graph.numComponents = numComponents;
results.graph.charPathLength = charPathLength;
results.graph.networkDiameter = networkDiameter;
results.graph.globalEfficiency = globalEfficiency;
results.graph.density = density;
results.graph.meanClusteringCoeff = meanClusteringCoeff;

results.stats.p_degree_early_g1_vs_g4 = p_degree_early;
results.stats.p_degree_late_g1_vs_g4 = p_degree_late;

if verbose
    fprintf('\nDegree centrality g1 vs g4, early CS trials p = %.4g\n', p_degree_early);
    fprintf('Degree centrality g1 vs g4, late CS trials p = %.4g\n', p_degree_late);
end

end
