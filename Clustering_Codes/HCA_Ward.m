function [Centroids,index,Cluster_Vectors,Cluster_n,Dendrogram_Total] = HCA_Ward(Distance_Metric,Cluster_range,Data)
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
Data1 = Data;
Data_Sum = Data;
w = 1; % Iteration count
Clusters_sum_num = num2cell(1:size(Data, 2)); % Cell array of initial clusters
Dendrogram_Total = zeros(size(Data, 2)-1, 3); % To store dendrogram info
Number_of_Clusters = 0; % Control variable for cluster number

% Main hierarchical clustering loop 
while Number_of_Clusters ~=1 %change this variable to change the number of clusters produced
% Calculate distances based on selected metric
SSE_Dist = cell(1,size(Data,2));%preallocation of variable
n = 1;%start n at 1

if contains("DTW", Distance_Metric)
    % Dynamic Time Warping (DTW) distance
    while n <= size(Data, 2)
        for e = 1:size(Data, 2)
            if ~isempty(Clusters_sum_num{n})
                for ii = 1:length(Clusters_sum_num{n})
                    SSE_Dist{n}(ii, e) = dtw(Data(~isnan(Data(:, n)),n), Data(~isnan(Data(:, e)),e));
                end
            else
                SSE_Dist{n}(ii, e) = 0;
            end
        end
        n = n + 1;
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
SSE_sqrt = cell(1, size(Data, 2));
for n = 1:size(Data, 2)
    for e = 1:size(Data, 2)
        SSE_sqrt{n}(1, e) = sqrt(2 * numel(Clusters_sum_num{n}) * numel(Clusters_sum_num{e}) / ...
                                  (numel(Clusters_sum_num{n}) + numel(Clusters_sum_num{e})));
    end
end

% Calculate the final SSE values for merging
SSE = zeros(size(Data, 2), size(Data, 2));
for n = 1:size(Data, 2)
    for e = 1:size(Data, 2)
        SSE(n, e) = SSE_sqrt{n}(1, e) * SSE_Dist{n}(1, e);
    end
end

% Rearrange SSE matrix and extract the minimum value (for merging clusters)
SSE_Rearranged = sort(SSE, 2);
SSE_Rearranged(:, 1) = []; % Remove first column with zeros
SSE_Min_Column = SSE_Rearranged(:, 1);

% Replace zeros with NaN
SSE_Min_Column(SSE_Min_Column == 0) = NaN;

% Find minimum SSE and its corresponding clusters
SSE_Min = min(SSE_Min_Column);
SSE_Min_Find = find(SSE_Min_Column == SSE_Min);

% Merge the two clusters identified by SSE_Min_Find
Clusters_sum_num{end + 1} = [Clusters_sum_num{SSE_Min_Find(1)}, Clusters_sum_num{SSE_Min_Find(2)}];

% Sum and average the merged clusters
Cluster_Sum = sum(Data(:, [SSE_Min_Find(1), SSE_Min_Find(2)]), 2);
Data_Sum(:, end + 1) = Cluster_Sum;
Data_Sum(:, [SSE_Min_Find(1), SSE_Min_Find(2)]) = 0;  % Remove merged columns

%This variable finds the time-series that match the columns from the
%SSE_Min_Find variable
Hierarchical_Vectors = Data1(:,Clusters_sum_num{end});

% Compute new cluster centroids as the mean of assigned data points
Cluster_Ave = zeros(size(Hierarchical_Vectors,1),1);
for i = 1:size(Hierarchical_Vectors,1)
    len = sum(~isnan(Hierarchical_Vectors(i,:)),'all');
    summation = nansum(Hierarchical_Vectors(i,:)); 
    Cluster_Ave(i,1) = summation/len;
end

Data(:,SSE_Min_Find) = -111111111111111111; 

%This adds the new averaged cluster to the end of the Matrix data
Data(:,end+1) = Cluster_Ave;

% Update Dendrogram with merge information
Dendrogram = [SSE_Min_Find(1),SSE_Min_Find(2), SSE_Min];
Dendrogram_Total(w, :) = Dendrogram;

% Mark the merged clusters as zero in the Clusters_sum_num list
Clusters_sum_num{SSE_Min_Find(1)} = 0;
Clusters_sum_num{SSE_Min_Find(2)} = 0;

% Update number of clusters
Number_of_Clusters = sum(cellfun(@(x) any(x ~= 0), Clusters_sum_num));

%Finds all the time-seris within each cluster
if Number_of_Clusters <= Cluster_range
    
    % Remove empty clusters and collect results
    Clusters_new = Clusters_sum_num(~cellfun(@(x) all(x == 0), Clusters_sum_num));
    
    Vectors = cell(1, Number_of_Clusters);
    for ii = 1:Number_of_Clusters
        Vectors{ii} = Data1(:, Clusters_new{ii});
    end

    % Create the index of clusters
    idx = zeros(size(Data1, 2), 1);
    for ii = 1:Number_of_Clusters
        idx(Clusters_new{ii}) = ii;
    end

    %Cluster_index{Number_of_Clusters,1} = idx;
    index{Number_of_Clusters,1} = Clusters_new;

    Cluster_Vectors{Number_of_Clusters,1} = Vectors;
    
    % Compute new cluster centroids as the mean of assigned data points
    Cluster_Ave = zeros(size(Data,1),Number_of_Clusters); %preallocate variable
    for ii = 1:Number_of_Clusters
        for i = 1:size(Cluster_Vectors{Number_of_Clusters,1}{ii},1)
            len = sum(~isnan(Cluster_Vectors{Number_of_Clusters,1}{ii}(i,:)),'all');
            summation = nansum(Cluster_Vectors{Number_of_Clusters,1}{ii}(i,:)); 
            ave(i,1) = summation/len;
        end
        Cluster_Ave(:,ii) = ave;
    end
    
    Centroids{Number_of_Clusters,1} = Cluster_Ave;
    
    % Number of participants in each cluster
    Cluster_n{Number_of_Clusters,1} = cellfun(@length, Clusters_new);

end

w = w+1;%w will increase by 1 after each run of the main while loop
end