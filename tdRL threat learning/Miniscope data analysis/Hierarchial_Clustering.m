function results = Hierarchial_Clustering(trace_shock_zs, clustnum, pre_baseline, maxNumClusters, writeClusterFiles, outputDir)
%% Hierarchical clustering of shock-aligned neural activity
% This function performs hierarchical clustering on shock-aligned z-scored
% neural activity traces and plots the mean activity of each cluster.
%
% @author:  Mingyang Wei, Shanghai Jiao Tong University, 2026
%
% Inputs:
%       - trace_shock_zs: z-scored neural activity aligned to shock
%                         neurons x timeBins x trials
%       - clustnum: number of clusters for final clustering
%       - pre_baseline: number of baseline frames before shock onset
%       - maxNumClusters: maximum number of clusters tested in elbow plot
%       - writeClusterFiles: whether to save cluster mean traces as files
%       - outputDir: output folder for cluster mean files
%
% Outputs:
%       - results: structure containing clustering outputs and cluster means
%             results.D
%             results.Z
%             results.dend
%             results.cophenetic_corr
%             results.wss_values
%             results.clustidx
%             results.savall
%             results.trace_shock_filtered_zs
%

%% Input defaults
if nargin < 2 || isempty(clustnum)
    clustnum = 4;
end

if nargin < 3 || isempty(pre_baseline)
    pre_baseline = 20;
end

if nargin < 4 || isempty(maxNumClusters)
    maxNumClusters = 10;
end

if nargin < 5 || isempty(writeClusterFiles)
    writeClusterFiles = false;
end

if nargin < 6 || isempty(outputDir)
    outputDir = pwd;
end

if writeClusterFiles && ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Filter shock-aligned z-scored traces
trace_shock_filtered_zs = [];
disp_results = 0;

for trail = 1:size(trace_shock_zs,3)
    for i = 1:12:size(trace_shock_zs,1)
        for j = 1:12
            try
                upsampled = interp(trace_shock_zs(i+j-1,:,trail),4);
                upsampled = smoothdata(upsampled,'gaussian','SmoothingFactor',0.85);
                downsampled = downsample(upsampled,4);
                trace_shock_filtered_zs(i+j-1,:,trail) = downsampled;

                if disp_results
                    subplot(3,4,j)
                    hold on
                    plot(1:size(trace_shock_zs,2),trace_shock_zs(i+j-1,:,trail))
                    plot(1:size(trace_shock_filtered_zs,2),trace_shock_filtered_zs(i+j-1,:,trail))
                    hold off
                end
            catch
                % out of bounds
            end
        end
        if disp_results
            pause(0.1)
            clf
        end
    end
end

%% Correlation averaged shock
D = squareform(pdist(mean(trace_shock_zs(:,1+pre_baseline:end,1:5),3),"correlation"));
Z = linkage(D,"ward");
figure;
dend = dendrogram(Z,0);
cophenetic_corr = cophenet(Z,D);

%% Elbow rule, see Bo Li 2023 for reference
wss_values = zeros(maxNumClusters, 1);

X = trace_shock_zs(:,1+pre_baseline:end,1:5);

for k = 1:maxNumClusters
    cluster_labels = cluster(Z, 'maxclust', k);
    % calculate WSS for this specific clustering
    wss_k = 0;
    for i = 1:k
        cluster_points = X(cluster_labels == i, :);
        cluster_center = mean(cluster_points, 1);
        wss_k = wss_k + sum(bsxfun(@minus, cluster_points, cluster_center).^2,"all");
    end
    wss_values(k) = wss_k;
end

% Plot WSS against number of clusters
figure;
plot(1:maxNumClusters, wss_values);
xlabel('Number of clusters');
ylabel('Total within-cluster sum of square');
title('Elbow method for determining number of clusters');

%% Cluster shock-aligned traces
savall = [];
clustidx = cluster(Z,maxclust=clustnum);  % first cluster index
color = ['k','r','y','g','b']; % since mod, 1-5 are rygbk: red yellow green blue black
meanall = []; % meanall/stdall: mean and std for each cluster in each trial
stdall = []; 
std_trace = std(trace_shock_zs,1,"all"); 
figure

for i = 1:clustnum
    subplot(2,4,i)
    hold on
    trace_use_all = mean(trace_shock_filtered_zs(:,:,1:5),3);
    plotmean = mean(trace_use_all(find(clustidx==i),:),1);
    plot(1:size(trace_shock_zs,2),plotmean,'Color','r')
    ylim([-1,3.5])
    savall = vertcat(savall,mean(trace_use_all(find(clustidx==i),:),1));
    title("num="+size(find(clustidx==i),1))
    hold off
    %xticks(pre_baseline)
    xticks([20,40])
    xticklabels({'US start','US end'})

    if writeClusterFiles
        writematrix(plotmean, fullfile(outputDir, 'Cluster'+string(i)+'.xls'));
    end
end

%% Save outputs
results = struct();
results.D = D;
results.Z = Z;
results.dend = dend;
results.cophenetic_corr = cophenetic_corr;
results.wss_values = wss_values;
results.clustidx = clustidx;
results.savall = savall;
results.trace_shock_filtered_zs = trace_shock_filtered_zs;
results.clustnum = clustnum;
results.pre_baseline = pre_baseline;
results.maxNumClusters = maxNumClusters;

end
