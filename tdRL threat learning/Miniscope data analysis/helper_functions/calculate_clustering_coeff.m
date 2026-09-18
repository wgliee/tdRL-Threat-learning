function clustering_coeffs = calculate_clustering_coeff(adjacency_matrix)
    n = size(adjacency_matrix, 1);  % 获取节点数
    clustering_coeffs = zeros(n, 1);  % 初始化聚类系数向量

    for i = 1:n
        neighbors = find(adjacency_matrix(i, :));  % 找到节点i的所有邻居
        if length(neighbors) > 1
            subgraph = adjacency_matrix(neighbors, neighbors);  % 邻居的邻接子矩阵
            possible_links = length(neighbors) * (length(neighbors) - 1) / 2;
            actual_links = sum(subgraph(:)) / 2;  % 子图中实际的边数，除以2是因为每条边被计算了两次
            clustering_coeffs(i) = 2 * actual_links / possible_links;
        else
            clustering_coeffs(i) = 0;  % 如果只有一个或没有邻居，则聚类系数为0
        end
    end
end
