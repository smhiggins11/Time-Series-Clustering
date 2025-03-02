function [Centroids,index,Vectors,Cluster_n,Initial_Centroids,empty_cluster,not_converged] = Kmeans(Distance_Metric,k,Data,Initial_Centroids)
% Created by Seth Higgins for PhD dissertation projects
% K-means clustering algorithm for time-series data
% 
% Inputs:
%   Distance_Metric   - Distance metric ('DTW' or 'Euclidean')
%   k                 - Number of clusters
%   Data              - Data matrix (time-series data as columns)
%   Initial_Centroids - (Optional) Predefined initial centroids
%
% Outputs:
%   Centroids         - Final cluster centroids
%   index             - Indices of data points assigned to each cluster
%   Vectors           - Clustered data points
%   Cluster_n         - Number of elements in each cluster
%   Initial_Centroids - Initial centroid indices
%   empty_cluster     - Number of times clusters were reinitialized due to empty clusters
%   not_converged     - Number of times clusters failed to converge within 50 iterations

%% Initialization
% If initial centroids are provided, use them; otherwise, generate random unique indices
if exist("Initial_Centroids",'var') == 1
    if isequal(size(Initial_Centroids), [size(Data,1), k])
        % Directly use the provided centroid matrix
        Centroids = Initial_Centroids;
        [~, Initial_Centroids] = ismember(Initial_Centroids', Data', "rows");
    elseif isequal(size(Initial_Centroids), [1, k])
        % Use provided indices to extract centroids from Data
        Centroids = Data(:, Initial_Centroids);
    else
        error("Initial_Centroids must be either a matrix of size [%d, %d] or a row vector of %d indices.", size(Data,1), k, k);
    end
else
    % Randomly select k unique indices and assign corresponding columns as centroids
    Initial_Centroids = randperm(size(Data,2), k);
    Centroids = Data(:, Initial_Centroids);
end

% Initialize counters for tracking empty clusters and non-converging clusters
empty_cluster = 0; 
not_converged = 0;

%% K-means Algorithm

Cluster_Converge_Ave = 1; % Convergence indicator (set to arbitrary nonzero value initially)
n = 0; % Iteration counter
while Cluster_Converge_Ave ~= 0 
    % Compute distances between each data point and cluster centroids  
    Cluster_Distance = zeros(k,size(Data,2));%preallocate variable
    
    if contains(Distance_Metric,"DTW")
        for ii = 1:k
            for a=1:size(Data,2)
                Cluster_Distance(ii,a) = dtw(Centroids(~isnan(Centroids(:,ii)),ii), Data(~isnan(Data(:,a)),a));
            end
        end
    elseif contains(Distance_Metric,"Euclidean")
        for ii = 1:k
            for a=1:size(Data,2)
                Cluster_Distance(ii,a) = norm(Centroids(:,ii) - Data(:,a))^2;
            end
        end
    end

    % Assign each data point to the nearest cluster centroid
    [~, Cluster_min_distance] = min(Cluster_Distance, [], 1);
    index = arrayfun(@(ii) find(Cluster_min_distance == ii), 1:k, 'UniformOutput', false);

    % Extract cluster data points
    Vectors = cell(1, k);
    for ii = 1:k
        Vectors{ii} = Data(:, index{ii});
    end

    % Compute new cluster centroids as the mean of assigned data points
    Cluster_Ave = cell(1,k); %preallocate variable
    for ii = 1:k
        for i = 1:size(Vectors{1,ii},1)
            len = sum(~isnan(Vectors{1,ii}(i,:)),'all');
            summation = nansum(Vectors{1,ii}(i,:)); 
            ave(i,1) = summation/len;
        end
        Cluster_Ave{ii} = ave;
    end

    % % Use DBA to compute the average for each cluster
    % for ii = 1:k
    %     average timeseries when they are different lengths
    %     if any(isnan(Vectors{1,ii}{1,1}),'all')
    %         if ~isempty(index{1,ii})
    %             Cluster_Ave{ii} = DBA(Vectors{1,ii}{1,1},Centroids(:,ii));
    %         else
    %             warning('Empty cluster detected, reinitializing centroid.');
    %             Initial_Centroids = randperm(size(Data,2));
    %             Initial_Centroids = Initial_Centroids(1:k);
    %             Centroids(:,ii) = Data(:,Initial_Centroids(ii));
    %         end
    %     end
    %     average timeseries when time-series are the same length \
    %     if all(~isnan(Vectors{1,ii}{1,1}),'all')
    %         Cluster_Ave{ii} = mean(Vectors{1,ii}{1,1},2);
    %     end
    % end

    % Check for convergence by comparing with previous cluster assignments
    if n > 0
        Cluster_Converge_Ave = mean(cellfun(@(x, y) dtw(x, y), index, Cluster_index_Previous));
    end

    % Update previous cluster assignments
    Cluster_index_Previous = index;

    % Update centroids for the next iteration
    for ii = 1:k
        Centroids(:, ii) = Cluster_Ave{ii};
    end

    % Check for duplicate centroids, which can cause instability
    for i = 1:k-1
        for j = i+1:k
            if isequal(Centroids(:, i), Centroids(:, j))
                error('Duplicate centroids detected between clusters %d and %d.', i, j);
            end
        end
    end

    n = n + 1; % Increment iteration counter
end

% Output number of elements in each cluster
Cluster_n = cellfun(@numel, index);

%sometimes clustering results leads to empty clusters. This will cause an
%error. below for loop will reset cluster centroids when there is an empty
%cluster
% Loop through each pair of time-series
for i = 1:k-1
    for j = i+1:k
        % Check if time-series i and j are identical
        if isequal(Centroids(:, i), Centroids(:, j))
            % Initial_Centroids = randperm(size(Data,2));
            % Initial_Centroids = Initial_Centroids(1:k);
            % Cluster_Centroids(:,ii) = Data(:,Initial_Centroids(ii));
            % n = 0;
            % fprintf('Time-series %d and %d are identical.\n', i, j);
            error('Time-series %d and %d are identical.\n', i, j);
            %break
        end
    end
end
dbstop if error

if n == 100 %if statment will tell me when the while loop has completed 50 cycle. At this point, the code may not converge and I may need to restart the code.

    Initial_Centroids = randperm(size(Data,2));
    Initial_Centroids = Initial_Centroids(1:k);
    uniqueCentroids = unique(Initial_Centroids);
    if length(uniqueCentroids) ~= length(Initial_Centroids)
        error('There are duplicate numbers in list')
    end
    
    Centroids = zeros(size(Data,1),k); 
    for ii = 1:k
        Centroids(:,ii) = Data(:,Initial_Centroids(ii));
    end
    n = 0;
    not_converged = not_converged+1;
    disp('Not converged')
end


end


