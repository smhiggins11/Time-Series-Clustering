function DI = dunn_index(k,Data, Cluster_Centroids, Cluster_index,Distance_Metric)
    % Function to compute the Dunn Index for evaluating clustering quality.
    % 
    % Inputs:
    %   k - Number of clusters
    %   Data - Time-series data (each column is a time series)
    %   Cluster_Centroids - Cluster centroids (each column is a centroid)
    %   Cluster_index - Cell array containing indices of data points in each cluster
    %   Distance_Metric - String specifying distance metric ('DTW' or 'Euclidean')
    %
    % Output:
    %   DI - Dunn Index (higher value indicates better clustering)
    
    % Initialize an index array to store cluster assignments
    idx = zeros(size(Data,2),1);
    % Assign each data point to its respective cluster based on
    % Cluster_index,,,,,,,,
    for i = 1:k
        for ii = 1:length(Cluster_index{1,i})
            idx(Cluster_index{1,i}(ii)) = i;
        end
    end
    
    % Precompute NaN removal for efficiency in DTW and Euclidean calculations
    cleaned_Centroids = cell(1, k);  % Cell array to store NaN-free centroids
    for ii = 1:k
        cleaned_Centroids{ii} = Cluster_Centroids(~isnan(Cluster_Centroids(:, ii)),ii);
    end


    % Initialize with a large value for minimum inter-cluster distance
    min_inter_distance = realmax;
    
    % Compute minimum distance between different cluster centroids
    if contains("DTW",Distance_Metric)
        % If using DTW, calculate the smallest DTW distance between cluster centroids
        for i = 1:k - 1
            for j = i+1:k
                
                distance = dtw(cleaned_Centroids{:,i}, cleaned_Centroids{:,j});
                % Update minimum inter-cluster distance if a smaller value is found
                if distance < min_inter_distance
                    min_inter_distance = distance;
                end
            end
        end
    elseif contains("Euclidean",Distance_Metric)
        % If using Euclidean distance, calculate the minimum squared distance between centroids
        for i = 1:k - 1
            for j = i+1:k
                distance = norm(Cluster_Centroids(:, i) - Cluster_Centroids(:, j))^2;
                if distance < min_inter_distance
                    min_inter_distance = distance;
                end
            end
        end
    end
    
    % Organize time-series data into clusters
    cluster_time_series = cell(1, k);
    for ii = 1:k
        cluster_time_series{ii} = Data(:, idx == ii);  % More efficient indexing
    end
    
    % Initialize cell array to store intra-cluster distances
    intra_distM = cell(k,1);
    %calculate maximum distance between any two time-series in the same cluster.
    if contains("DTW",Distance_Metric)

        cleaned_vectors = cell(1,k);
        for ii = 1:k
            for i = 1:size(cluster_time_series{ii},2)
                cleaned_vectors{ii}{i} = cluster_time_series{ii}(~isnan(cluster_time_series{1,ii}(:,i)),i);
            end
        end

        % Compute DTW distances between all time-series pairs in each cluster
        for ii = 1:k
            for i = 1:size(cluster_time_series{1,ii},2)
                for a = 1:size(cluster_time_series{1,ii},2)
                    intra_distM{ii,1}(a,i) = dtw(cleaned_vectors{ii}{i},cleaned_vectors{ii}{a});
                end
            end
        end
    elseif contains("Euclidean",Distance_Metric)
        % Compute Euclidean distances between all time-series pairs in each cluster
        for ii = 1:k
            for i = 1:size(cluster_time_series{1,ii},2)
                for a = 1:size(cluster_time_series{1,ii},2)
                    intra_distM{ii,1}(a,i) = norm(cluster_time_series{1,ii}(:,i) - cluster_time_series{1,ii}(:,a))^2;
                end
            end
        end
    end
      
    % Directly extract the maximum intra-cluster distance from each cluster and get the largest across all clusters
    max_intra_distance = max(cellfun(@(x) max(x(:)), intra_distM));

    % Compute the Dunn Index as the ratio of min inter-cluster distance to max intra-cluster distance
    DI = min_inter_distance / max_intra_distance;
end