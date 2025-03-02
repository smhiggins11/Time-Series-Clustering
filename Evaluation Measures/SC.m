%% Silhouette Coefficient
function Silhouette_Score = SC(k,Cluster_Vectors,Cluster_Centroids,Cluster_n,Data,Distance_Metric)

% Function to compute the Silhouette Coefficient for clustering analysis
% Inputs:
%   k - Number of clusters
%   Cluster_Vectors - Cell array containing cluster data
%   Cluster_Centroids - Matrix of cluster centroids
%   Cluster_n - Number of elements in each cluster
%   Data - Original dataset
%   Distance_Metric - String specifying distance metric ('DTW' or 'Euclidean')
% Output:
%   Silhouette_Score - Computed silhouette score for the clustering

%% Explanation of the calculations
% calculate a = average distance of i to the points in the same cluster
% calculate b = min(average distance of i to points in another cluster)
% Silhouette coefficient = 1 - a/b 
% if a < b, (or s = b/a - 1 if a > or equal to b, not the usual case)
% typically between 0 and 1
% the closer to 1 the better

%This calculates (a) variable in the silhouette equation.
%This calculates the distance of each vector to other vectors within the
%same cluster.

% Compute intra-cluster distances (a) - Average distance of a point to other points in the same cluster
Cluster_Distance_Within = struct([]);
if contains("DTW",Distance_Metric)
    for ii = 1:k
        n = 1;
        while n<size(Cluster_Vectors{1,ii},2)+1
            for a = 1:size(Cluster_Vectors{1,ii},2) 
                vector_a = Cluster_Vectors{1,ii}(:,n);
                vector_b = Cluster_Vectors{1,ii}(:,a);
                % Compute DTW distance between vectors within the same cluster
                Cluster_Distance_Within(ii).Silhouette(a,n) = dtw(vector_a(~isnan(vector_a)), vector_b(~isnan(vector_b)));
            end
            n = n+1;
        end
    end
end

if contains("Euclidean",Distance_Metric)
    for ii = 1:k
        n = 1;
        while n<size(Cluster_Vectors{1,ii},2)+1
            for a = 1:size(Cluster_Vectors{1,ii},2) 
                vector_a = Cluster_Vectors{1,ii}(:,n);
                vector_b = Cluster_Vectors{1,ii}(:,a);
                % Compute squared Euclidean distance between vectors in the same cluster
                Cluster_Distance_Within(ii).Silhouette(a,n) = norm(vector_a(~isnan(vector_a)) - vector_b(~isnan(vector_b)))^2;
            end
            n = n+1;
        end
    end
end

% Sort distances and remove self-distance (which is zero)
for ii = 1:k
    for a = 1:size(Cluster_Distance_Within(ii).Silhouette,2)
        Cluster_Distance_Within(ii).Silhouette(a,:) = sort(Cluster_Distance_Within(ii).Silhouette(a,:));
    end
end

% Remove zero distances (self-distance)
for ii = 1:k
    if length(Cluster_Distance_Within(ii).Silhouette(:,1)) > 1
        Cluster_Distance_Within(ii).Silhouette(:,1) = [];
    elseif length(Cluster_Distance_Within(ii).Silhouette(:,1)) == 1
        Cluster_Distance_Within(ii).Silhouette(:,1) = Cluster_Distance_Within(ii).Silhouette(:,1);
    end
end

% Compute average intra-cluster distance (a)
for ii = 1:k
    for a = 1:length(Cluster_Distance_Within(ii).Silhouette)
        Cluster_Distance_Within(ii).Silhouette_Ave(a,:) = nanmean(Cluster_Distance_Within(ii).Silhouette(a,:));
    end
end

% Compute inter-cluster distances (b) - Distance of each point to other cluster centroids
Cluster_Distance = struct([]);
if contains("DTW",Distance_Metric)
    for ii = 1:k
        for b = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
                vector = Cluster_Vectors{1,ii}(:,a);
                centroid = Cluster_Centroids(:,b);
                % Compute DTW distance to other cluster centroids
                Cluster_Distance(ii).Between(b,a) = dtw(vector(~isnan(vector)), centroid(~isnan(centroid)));
            end
        end
    end
end

if contains("Euclidean",Distance_Metric)
    for ii = 1:k
        for b = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
                vector = Cluster_Vectors{1,ii}(:,a);
                centroid = Cluster_Centroids(:,b);
                % Compute Euclidean distance to other cluster centroids
                Cluster_Distance(ii).Between(b,a) = norm(vector(~isnan(vector)) - centroid(~isnan(centroid)))^2;
            end
        end
    end
end

% Remove distances to the same cluster (not needed)
for ii = 1:k
    Cluster_Distance(ii).Between(ii,:) = [];
end

% Compute minimum distance to another cluster (b)
if k > 2
    for ii = 1:k
        Cluster_Distance(ii).Min = min(Cluster_Distance(ii).Between);
    end
end

% Compute silhouette score for each point
Silhouette = zeros(size(Data,2),1);
b = 1;
for ii = 1:k
    for a = 1:length(Cluster_Distance_Within(ii).Silhouette_Ave)
        a_dist = Cluster_Distance_Within(ii).Silhouette_Ave(a); % Intra-cluster distance (a)
        b_dist = Cluster_Distance(ii).Min(a); % Nearest other cluster distance (b)
        % Compute silhouette coefficient
        Silhouette(b) = (b_dist - a_dist) / max(a_dist, b_dist);
        b = b+1;
    end
end

Silhouette_Score = mean(Silhouette);

end