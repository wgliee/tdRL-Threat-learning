function [] = rowdata_columnumber_plot_omitnan(inputmat,toscatter)
% input: [—————
%         ——x——
%         —————]    num*dim

if nargin==1
    toscatter=0;
end

hold on
rowsize=size(inputmat,1);
colsize=size(inputmat,2);

if toscatter
    for i=1:rowsize
    scatter(1:colsize,inputmat(i,:))
    end
end

% errorbar
% errorbar
standard_deviation = nanstd(inputmat, 0, 1); % Compute the standard deviation along the 1st dimension while ignoring NaNs.
mean_values = nanmean(inputmat, 1); % Compute the mean along the 1st dimension while ignoring NaNs.

errorbar(mean_values,standard_deviation./sqrt(rowsize))

hold off
end