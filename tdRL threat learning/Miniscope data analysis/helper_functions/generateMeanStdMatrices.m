function [meanMatrix, stdMatrix] = generateMeanStdMatrices(tone_dataset)
    % tone_dataset： dims*n*tones
    % Source: Neuron 2024 Methods
    %   Excitability mediates allocation of pre-configured
    %   ensembles to a hippocampal engram supporting
    %   contextual conditioned threat in mice
    % by wmy in 18-06-2024

    % 获取tones的数量
    numTones = size(tone_dataset, 3);
    % 初始化mean和std矩阵
    dims = size(tone_dataset, 1);
    meanMatrix = zeros(dims, dims, numTones);
    stdMatrix = zeros(dims, dims, numTones);
    
    % 对每个tone计算mean和std矩阵
    for tones = 1:numTones
        temp = tone_dataset(:,:,tones)';
        dims = size(temp, 2);
        allConnectivityMatrices = zeros(dims, dims, 1000); % 初始化存储1000个相关矩阵的容器
        
        % 进行1000次随机打乱和相关性计算
        for iteration = 1:1000
            shuffledTemp = temp;
            for col = 1:dims
                % 随机选择位移量k，范围从0到numRows-1
                k = randi(dims) - 1;  % MATLAB的randi生成1到numRows的随机数，因此需要减1
                % 对该列进行循环位移
                shuffledTemp(:, col) = circshift(temp(:, col), k);
            end
            % 计算打乱后的相关性矩阵
            connectivityMat = corr(shuffledTemp, shuffledTemp);
            % 储存每次生成的相关性矩阵
            allConnectivityMatrices(:,:,iteration) = connectivityMat;
        end
        
        % 计算所有相关性矩阵的平均值和标准差
        meanMatrix(:,:,tones) = mean(allConnectivityMatrices, 3);
        stdMatrix(:,:,tones) = std(allConnectivityMatrices, 0, 3);
    end
end
