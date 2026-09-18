function [idx_dat1, idx_dat2] = random_subsample_getidx(input_dat1, input_dat2)
% RANDOM_SUBSAMPLE Randomly subsamples the larger dataset to match the size of the smaller one.
%   Both input datasets are expected to have dimensions: (downsample dim x same dim)
%   Returns indices for the subsampled data from the original datasets.

% Determine which dataset is larger
if size(input_dat1, 1) < size(input_dat2, 1)
    larger_dat = input_dat2;
    smaller_dat = input_dat1;
    switchFlag = true;
else
    larger_dat = input_dat1;
    smaller_dat = input_dat2;
    switchFlag = false;
end

% Randomly subsample the larger dataset
numFrames = size(smaller_dat, 1);
randomIdx = randperm(size(larger_dat, 1), numFrames);

% Assign outputs based on the original order of the inputs
if switchFlag
    idx_dat1 = 1:size(smaller_dat, 1); % since smaller_dat corresponds to input_dat1
    idx_dat2 = randomIdx;
else
    idx_dat1 = randomIdx;
    idx_dat2 = 1:size(smaller_dat, 1); % since smaller_dat corresponds to input_dat2
end



end