function CHIndex = CH_index(Matrix,cluster_index,Distance_Metric)
% CH_index: Computes the Calinski-Harabasz Index for clustering validation.
%
% INPUTS:
%   - Matrix: A (m x n) matrix where each column is a data point (time series).
%   - cluster_index: Either a vector (1 x n) of cluster assignments or a cell array 
%                    where each cell contains indices of data points belonging to a cluster.
%   - Distance_Metric: String specifying the distance metric ('Euclidean' or 'DTW').
%
% OUTPUT:
%   - CHIndex: The Calinski-Harabasz index, used to evaluate clustering quality.
%
% The CH index is computed as:
%       CH = (BSS / (k - 1)) / (WSS / (n - k))
%   where:
%       - BSS (Between-cluster sum of squares) measures separation between clusters.
%       - WSS (Within-cluster sum of squares) measures compactness within clusters.
%
% This function supports both Euclidean and Dynamic Time Warping (DTW) distance metrics.

%% Convert cluster_index to an index vector (idx)
if isa(cluster_index, 'double')
    % Directly assign if already a numerical vector
    idx = cluster_index;
elseif isa(cluster_index, 'cell')
    % Convert cell array to a numerical index vector
    idx = zeros(1, size(Matrix, 2));
    for ii = 1:length(cluster_index)
        idx(cell2mat(cluster_index(ii))) = ii;
    end
end

% Number of clusters (k) and total data points (n)
k = max(idx); % Number of clusters
n = size(Matrix, 2); % Number of data points

% Compute the overall mean of all data points (ignoring NaNs)
overallMean = nanmean(Matrix,2); % needs to be 1X(length of vector)

% Initialize Within-Cluster Sum of Squares (WSS) and Between-Cluster Sum of Squares (BSS)
WSS = 0;
BSS = 0;
%% Compute CH Index using Euclidean distance
if contains(Distance_Metric,'Euclidean')
    for i = 1:k
        % Points in cluster i
        clusterPoints = Matrix(:,idx == i);
        clusterMean = mean(clusterPoints,2); % needs to be 1X(length of vector)
        
        % Within-cluster sum of squares
        WSS = WSS + sum(sum((clusterPoints - clusterMean).^2));
        
        % Between-cluster sum of squares
        BSS = BSS + size(clusterPoints, 2) * sum((clusterMean - overallMean).^2);
    end
%% Compute CH Index using Dynamic Time Warping (DTW) distance
elseif contains('DTW',Distance_Metric)
    for i = 1:k
        % Points in cluster i
        clusterPoints = Matrix(:,idx == i);
        
        % Calculate cluster mean (using DTW)
        clusterMean = nanmean(clusterPoints,2); % Euclidean mean, could be replaced with a DTW-based centroid if needed
        
        % Within-cluster sum of squares (WSS) using DTW distance
        for j = 1:size(clusterPoints, 2)
            cluster = clusterPoints(:,j);
            WSS = WSS + dtw(cluster(~isnan(cluster)), clusterMean(~isnan(clusterMean))); % Replace with DTW distance
        end
        
        % Between-cluster sum of squares (BSS) using DTW distance
        for j = 1:size(clusterPoints, 2)
            BSS = BSS + size(clusterPoints, 2) *  dtw(clusterMean(~isnan(clusterMean)), overallMean(~isnan(overallMean))); % Replace with DTW distance
        end
    end
end

% Calinski-Harabasz Index
CHIndex = (BSS / (k - 1)) / (WSS / (n - k));