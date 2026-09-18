function [newtrace] = downsample_trace(inputmat, downsample_times)
% downsample_trace downsample the 2nd dimension of the input matrix.
% inputmat: a 2-D or 3-D matrix with dimensions [neurons, time, slices]
% downsample_times: an integer factor by which to downsample
% Output is a matrix with the 2nd dimension downsampled.

% Get the dimensions of the input matrix.
[neurons, frames, slices] = size(inputmat);

% Check if the 2nd dimension can be downsampled by the given factor.
assert(mod(frames, downsample_times) == 0, "Second dimension cannot be evenly divided by downsample_times!");

new_frames = frames / downsample_times;

% Initialize output
newtrace = zeros(neurons, new_frames, slices);

for k = 1:slices
    for i = 1:new_frames
        start_idx = 1 + (i-1)*downsample_times;
        end_idx = i*downsample_times;

        downsampled_temp = mean(inputmat(:, start_idx:end_idx, k), 2);
        newtrace(:, i, k) = downsampled_temp;
    end
end

end