function [Centroids,index,Cluster_Vectors,Cluster_n,Dendrogram_Total] = HCA_Ward_Faster(Distance_Metric,Cluster_range,Data)
% HCA_Ward performs Hierarchical Clustering using the Ward method.
% 
% Inputs:
%   Distance_Metric: Metric to calculate distances ('DTW' or 'Euclidean')
%   Cluster_range: Desired range of clusters (stopping criteria)
%   Data: Matrix of time series data (columns are subjects)
%
% Outputs:
%   Centroids: The centroids (mean) for each cluster
%   index: Index of subjects assigned to each cluster
%   Cluster_Vectors: Vectors for each cluster
%   Cluster_n: The number of elements in each cluster
%   Dendrogram_Total: The full hierarchical clustering results

% Initialize variables
Original_Data = Data;
%Data_Sum = Data;
w = 1; % Iteration count
Clusters_sum_num = num2cell(1:size(Data, 2)); % Cell array of initial clusters
Clusters_sum_num(end+1:end+size(Data, 2)-1) = {[]};
Dendrogram_Total = zeros(size(Data, 2)-1, 3); % To store dendrogram info
Number_of_Clusters = 0; % Control variable for cluster number
index = cell(Cluster_range,1);
Cluster_Vectors = cell(Cluster_range,1);
Centroids = cell(Cluster_range,1);
Cluster_n = cell(Cluster_range,1);

% Main hierarchical clustering loop 
while Number_of_Clusters ~=1 %change this variable to change the number of clusters produced
% Calculate distances based on selected metric
n = 1;%start n at 1
if contains("DTW", Distance_Metric)
    N = size(Data, 2);
    SSE_Dist = zeros(N, N);
    for n = 1:N
        if ~isempty(Clusters_sum_num{n})
            idx_n = n;
            x = Data(~isnan(Data(:, idx_n)), idx_n);
            for e = 1:N
                y = Data(~isnan(Data(:, e)), e);
                SSE_Dist(n, e) = dtw(x, y);
            end
        end
    end
elseif contains("Euclidean", Distance_Metric)
    % Euclidean distance
    while n <= size(Data, 2)
        for e = 1:size(Data, 2)
            if ~isempty(Clusters_sum_num{n})
                for ii = 1:length(Clusters_sum_num{n})
                    SSE_Dist{n}(ii, e) = norm(Data(:, n) - Data(:, e));
                end
            else
                SSE_Dist{n}(ii, e) = 0;
            end
        end
        n = n + 1;
    end
end

% Calculate the adjusted SSE for merging clusters
N = size(Data, 2);
cluster_sizes = zeros(1, N);
for i = 1:N
    cluster_sizes(i) = numel(Clusters_sum_num{i});
end

% Use meshgrid to create a matrix of combinations
[A, B] = meshgrid(cluster_sizes, cluster_sizes);
SSE_sqrt = sqrt(2 .* A .* B ./ (A + B + eps)); % +eps avoids divide-by-zero

% Calculate the final SSE values for merging
SSE = SSE_sqrt .* SSE_Dist;

% Rearrange SSE matrix and extract the minimum value (for merging clusters)
SSE_Rearranged = sort(SSE, 2);
SSE_Rearranged(:, 1) = []; % Remove first column with zeros
SSE_Min_Column = SSE_Rearranged(:, 1);

% Replace zeros with NaN
SSE_Min_Column(SSE_Min_Column == 0) = NaN;

% Find minimum SSE and its corresponding clusters
SSE_Min = min(SSE_Min_Column);
SSE_Min_Find = find(SSE_Min_Column == SSE_Min);

new_cluster = [Clusters_sum_num{SSE_Min_Find(1)}, Clusters_sum_num{SSE_Min_Find(2)}];
Clusters_sum_num{size(Data,2)+1} = new_cluster;

% Compute the average time series of the merged cluster (ignoring NaNs)
merged_data = Original_Data(:, new_cluster);
Cluster_Ave = nanmean(merged_data, 2); % This gives you the centroid

% Sum and average the merged clusters
% Cluster_Sum = sum(Data(:, [SSE_Min_Find(1), SSE_Min_Find(2)]), 2);
% Data_Sum(:, end + 1) = Cluster_Sum;
% Data_Sum(:, [SSE_Min_Find(1), SSE_Min_Find(2)]) = 0;  % Remove merged columns

Data(:,SSE_Min_Find) = NaN; 

% Append the averaged cluster to Data
Data(:, end + 1) = Cluster_Ave;

% Update Dendrogram with merge information
Dendrogram = [SSE_Min_Find(1),SSE_Min_Find(2), SSE_Min];
Dendrogram_Total(w, :) = Dendrogram;

% Mark the merged clusters as zero in the Clusters_sum_num list
Clusters_sum_num{SSE_Min_Find(1)} = {};
Clusters_sum_num{SSE_Min_Find(2)} = {};

% Update number of clusters
Number_of_Clusters = sum(cellfun(@(x) ~isempty(x), Clusters_sum_num));

%Finds all the time-seris within each cluster
if Number_of_Clusters <= Cluster_range
    
    % Remove empty clusters and collect results
    Clusters_new = Clusters_sum_num(~cellfun(@isempty, Clusters_sum_num));
    
    Vectors = cell(1, Number_of_Clusters);
    for ii = 1:Number_of_Clusters
        Vectors{ii} = Original_Data(:, Clusters_new{ii});
    end

    % Create the index of clusters
    idx = zeros(size(Original_Data, 2), 1);
    for ii = 1:Number_of_Clusters
        idx(Clusters_new{ii}) = ii;
    end

    %Cluster_index{Number_of_Clusters,1} = idx;
    index{Number_of_Clusters,1} = Clusters_new;

    Cluster_Vectors{Number_of_Clusters,1} = Vectors;
    
    % Compute new cluster centroids as the mean of assigned data points
    Cluster_Ave = zeros(size(Data,1),Number_of_Clusters); %preallocate variable
    for ii = 1:Number_of_Clusters
        Cluster_Ave(:,ii) = nanmean(Cluster_Vectors{Number_of_Clusters,1}{1,ii},2); % This gives you the centroid
    end
    
    Centroids{Number_of_Clusters,1} = Cluster_Ave;
    
    % Number of participants in each cluster
    Cluster_n{Number_of_Clusters,1} = cellfun(@length, Clusters_new);

end

w = w+1;%w will increase by 1 after each run of the main while loop
end
end