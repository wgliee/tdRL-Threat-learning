function output = timeSeriesAverage(inputSeries, n)
    % 获取输入序列的长度
    t = length(inputSeries);
    % 初始化输出序列
    output = zeros(1, t);
    
    % 对于输入序列的每个元素
    for i = 1:t
        % 计算平均值的开始和结束索引
        startIdx = max(1, i - n);
        endIdx = min(t, i + n-1);
        
        % 计算平均值
        output(i) = mean(inputSeries(startIdx:endIdx));
    end
end
