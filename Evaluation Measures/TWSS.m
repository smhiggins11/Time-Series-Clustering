function Within_Cluster_Total_Sum = TWSS(k,Cluster_Centroids,Cluster_Vectors,Distance_Metric)
%% Total Within Sum of Squares
% This function computes the Total Within Cluster Sum of Squares (TWSS) 
% based on the given distance metric (DTW or Euclidean) for evaluating 
% clustering performance.

% Inputs:
%   k - Number of clusters
%   Cluster_Centroids - Matrix containing the centroids of each cluster (one per column)
%   Cluster_Vectors - Cell array containing the time-series data for each cluster (in the same structure as the centroids)
%   Distance_Metric - String specifying the distance metric ('DTW' or 'Euclidean')

% Output:
%   Within_Cluster_Total_Sum - The total within-cluster sum of squares for the clustering

%% Initialize structure to hold distances for each cluster
Within_Cluster = struct('Distance',cell(1,k));%preallocate variable

% If using DTW distance metric
if contains("DTW",Distance_Metric)

    % Precompute NaN removal for efficiency in DTW and Euclidean calculations
    cleaned_Centroids = cell(1, k);  % Cell array to store NaN-free centroids
    for ii = 1:k
        cleaned_Centroids{ii} = Cluster_Centroids(~isnan(Cluster_Centroids(:, ii)),ii);
    end
    
    cleaned_vectors = cell(1,k);
    for ii = 1:k
        for i = 1:size(Cluster_Vectors{ii},2)
            cleaned_vectors{ii}{i} = Cluster_Vectors{ii}(~isnan(Cluster_Vectors{1,ii}(:,i)),i);
        end
    end

    % Loop through each cluster
    for ii = 1:k
        % Loop through each time-series vector in the current cluster
        for a = 1:size(cleaned_vectors{ii},2)
            % Compute DTW distance between the vector and the centroid 
            Within_Cluster(ii).Distance(a,1) = dtw(cleaned_vectors{ii}{a}, cleaned_Centroids{ii});
        end
    end
% If using Euclidean distance metric
elseif contains("Euclidean",Distance_Metric)
    % Loop through each cluster
    for ii = 1:k
        % Loop through each time-series vector in the current cluster
        for a = 1:size(Cluster_Vectors{1,ii},2)
            % Compute Euclidean squared distance between the vector and the centroid
            Within_Cluster(ii).Distance(a,1) = norm(Cluster_Vectors{1,ii}(:,a) - Cluster_Centroids(:,ii))^2;
        end
    end
end

%% Calculate the sum of distances within each cluster
Within_Cluster_Sum = zeros(1,k); % Preallocate array to store the sum of distances for each cluster
for ii = 1:k
    % Sum the distances for each vector in the current cluster
    Within_Cluster_Sum(ii) = sum(Within_Cluster(ii).Distance);
end

%% Compute the total within-cluster sum of squares
Within_Cluster_Total_Sum = sum(Within_Cluster_Sum); % Sum of all clusters' sums of squared distances
end