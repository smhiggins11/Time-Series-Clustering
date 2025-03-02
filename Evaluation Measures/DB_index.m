function DBI = DB_index(Matrix, centroids, cluster_index, Distance_Metric)
% DB_index computes the Davies-Bouldin Index (DBI) for clustering validation.
%
% Inputs:
%   - Matrix: Data matrix (features x samples)
%   - centroids: Cluster centroids (features x clusters)
%   - cluster_index: Cluster assignments (either as a double array or cell array)
%   - Distance_Metric: String specifying the distance metric ('Euclidean' or 'DTW')
%
% Output:
%   - DBI: Davies-Bouldin Index (lower values indicate better clustering)

% Convert cluster_index to numerical form if needed
if isa(cluster_index, 'double')
    idx = cluster_index;
elseif isa(cluster_index, 'cell')
    idx = zeros(1, size(Matrix, 2));
    for ii = 1:length(cluster_index)
        idx(cluster_index{ii}) = ii; % Assign cluster indices
    end
else
    error('cluster_index must be a double array or a cell array');
end

% Number of clusters
k = max(idx);

% Initialize intra-cluster (S) and inter-cluster (D) distances
S = zeros(k, 1); 
D = inf(k, k); % Set diagonal to inf to avoid division by zero

% Compute intra-cluster distances (S)
if contains(Distance_Metric, 'Euclidean')
    for i = 1:k
        clusterPoints = Matrix(:, idx == i);
        clusterCentroid = centroids(:, i);
        S(i) = mean(vecnorm(clusterPoints - clusterCentroid, 2, 1)); % Faster Euclidean norm
    end
elseif contains(Distance_Metric, 'DTW')
    for i = 1:k
        clusterPoints = Matrix(:, idx == i);
        clusterCentroid = centroids(:, i);
        dtwDistances = arrayfun(@(p) dtw(clusterPoints(~isnan(clusterPoints(:, p)),p), clusterCentroid(~isnan(clusterCentroid))), 1:size(clusterPoints, 2));
        S(i) = nanmean(dtwDistances); % Handle NaNs in DTW calculation
    end
else
    error('Unsupported distance metric. Use "Euclidean" or "DTW".');
end

% Compute inter-cluster distances (D)
if contains(Distance_Metric, 'Euclidean')
    for i = 1:k-1
        for j = i+1:k
            D(i, j) = norm(centroids(:, i) - centroids(:, j));
            D(j, i) = D(i, j); % Symmetric matrix
        end
    end
elseif contains(Distance_Metric, 'DTW')
    for i = 1:k-1
        for j = i+1:k
            D(i, j) = dtw(centroids(~isnan(centroids(:, i)),i), centroids(~isnan(centroids(:, j)),j));
            D(j, i) = D(i, j);
        end
    end
end

% Compute R_ij = (S_i + S_j) / D_ij
R = (S + S') ./ D;
R(isinf(R)) = 0; % Handle division by zero

% Compute Davies-Bouldin Index (DBI)
DBI = mean(max(R, [], 2)); % Max over columns, mean over rows

end