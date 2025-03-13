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


if k > 1 % Cannot calculate Silhouette Coefficient for one cluster. Number of clusters needs to be greater than 1
    if contains(Distance_Metric,"DTW")
        Cluster_Distance_Within = struct([]);
        for ii = 1:k
            n = 1;
            while n<size(Cluster_Vectors{1,ii},2)+1
                for a = 1:size(Cluster_Vectors{1,ii},2)
                    Cluster_Distance_Within(ii).Silhouette(a,n) = dtw(Cluster_Vectors{1,ii}(:,n), Cluster_Vectors{1,ii}(:,a));
                end
                n = n+1;
            end
        end
    elseif contains(Distance_Metric,"Euclidean")
        Cluster_Distance_Within = struct([]);
        for ii = 1:k
            n = 1;
            while n<size(Cluster_Vectors{1,ii},2)+1
                for a = 1:size(Cluster_Vectors{1,ii},2)
                    Cluster_Distance_Within(ii).Silhouette(a,n) = norm(Cluster_Vectors{1,ii}(:,n) - Cluster_Vectors{1,ii}(:,a))^2;
                end
                n = n+1;
            end
        end
    end

    %This sorts the values in each row from smallest to greatest number.
    %This is so I can delete the column with all the zeros. 
    %Zeros represent the distance between the same vector, which we do not
    %need.
    for ii = 1:k
        for a = 1:size(Cluster_Distance_Within(ii).Silhouette,2)
            Cluster_Distance_Within(ii).Silhouette(a,:) = sort(Cluster_Distance_Within(ii).Silhouette(a,:));
        end
    end
    
    %This deletes the column with the zeros
    for ii = 1:k
        if length(Cluster_Distance_Within(ii).Silhouette(:,1)) > 1
            Cluster_Distance_Within(ii).Silhouette(:,1) = [];
        elseif length(Cluster_Distance_Within(ii).Silhouette(:,1)) == 1
            Cluster_Distance_Within(ii).Silhouette(:,1) = Cluster_Distance_Within(ii).Silhouette(:,1);
        end
    end
    
    %This calculates the average of for each row vector
    for ii = 1:k
        for a = 1:length(Cluster_Distance_Within(ii).Silhouette)
            Cluster_Distance_Within(ii).Silhouette_Ave(a,:) = mean(Cluster_Distance_Within(ii).Silhouette(a,:));
        end
    end
    
    %This calculates the distance between each vector within each cluster to
    %the centroid of the other clusters. This will help determine the nearest
    %centroid each vector is closest two (other than the centroid they are already in)
    if contains(Distance_Metric,"DTW")
        Cluster_Distance = struct([]);
        for ii = 1:k
            for b = 1:k
                for a = 1:size(Cluster_Vectors{1,ii},2)
                    Cluster_Distance(ii).Between(b,a) = dtw(Cluster_Vectors{1,ii}(:,a), Cluster_Centroids(:,b));
                end
            end
        end
    elseif contains(Distance_Metric,"Euclidean")
        Cluster_Distance = struct([]);
        for ii = 1:k
            for b = 1:k
                for a = 1:size(Cluster_Vectors{1,ii},2)
                    Cluster_Distance(ii).Between(b,a) = norm(Cluster_Vectors{1,ii}(:,a) - Cluster_Centroids(:,b))^2;
                end
            end
        end
    end
    %The above section also calculated the distance between each vector and
    %their own centroid which we do not need. This section will delete those
    %rows because we do not need them. 
    for ii = 1:k
        Cluster_Distance(ii).Between(ii,:) = [];
    end
    
    %This calculates the minimum distance values within each vector distance
    %measures.
    if k > 2
        for ii = 1:k
            Cluster_Distance(ii).Min = min(Cluster_Distance(ii).Between);
        end
    end
    
    %This section subtracts the minimum values just calculated with the
    %Cluster_Dist.Between variable. A zero represents a vector that is closest
    %to that centroid.
    if k > 2
        for ii = 1:k
            for a = 1:size(Cluster_Distance(ii).Between,1)
                Cluster_Distance(ii).Sim(a,:) = Cluster_Distance(ii).Between(a,:) - Cluster_Distance(ii).Min;
            end
        end
    end
    
    %In order to have this code calculate silhouette coefficients for any
    %number of centroid k, I created the Clus variable, which creates a series
    %of number based on how many centroid there are.
    Cluster_k = zeros(k,k);
    for ii = 1:k
        Cluster_k(ii,1:k) = 1:k;
    end
    
    %This will eliminate the centroid number of that centroid. 
    %So with in the first column which is centroid 1, the 1 will be changed to
    %a zero and centroid two the two will be changed to a zero and so on.if k > 2
    for ii = 1:k
        for a = 1:k
            if Cluster_k(ii,a) == ii
                Cluster_k(ii,a) = 0;
            end
        end
    end
    
    %I want to delete the zeros within this variable so I start by sorting the
    %rows
    for ii = 1:k
        Cluster_k(ii,:) = sort(Cluster_k(ii,:));
    end
    
    
    %This deletes the column with the zeros
    Cluster_k(:,1) = [];
    
    if k > 2
        for ii = 1:k
            for a = 1:size(Cluster_Distance(ii).Sim,1) 
                for b = 1:size(Cluster_Distance(ii).Sim(a,:),2)
                    if Cluster_Distance(ii).Sim(a,b) == 0
                        Cluster_Distance(ii).Sim(a,b) = Cluster_k(ii,1);
                    elseif Cluster_Distance(ii).Sim(a,b) > 0
                        Cluster_Distance(ii).Sim(a,b) = Cluster_k(ii,2);
                    end
                end
            end
        end
    end
    
    %This section calculates the distance between each individual vector in one
    %cluster and to each vector in the other clusters.
    Silhouette = struct([]);
    for ii = 1:k
        for b = 1:k
            Silhouette(b,ii).Cluster = zeros(Cluster_n(1,ii),Cluster_n(1,b));
        end
    end
    
    if contains(Distance_Metric,"DTW")
        for ii = 1:k
            n = 1;
            while n<size(Cluster_Vectors{1,ii},2)+1
                for b = 1:k
                    for a = 1:size(Cluster_Vectors{1,b},2)
                        Silhouette(b,ii).Cluster(n,a) = dtw(Cluster_Vectors{1,ii}(:,n), Cluster_Vectors{1,b}(:,a));
                    end
                end
                n = n+1;
            end
        end
    elseif contains(Distance_Metric,"Euclidean")
        for ii = 1:k
            n = 1;
            while n<size(Cluster_Vectors{1,ii},2)+1
                for b = 1:k
                    for a = 1:size(Cluster_Vectors{1,b},2)
                        Silhouette(b,ii).Cluster(n,a) = norm(Cluster_Vectors{1,ii}(:,n) - Cluster_Vectors{1,b}(:,a))^2;
                    end
                end
                n = n+1;
            end
        end
    end

    
    %Since I only care about the distances in vectors to the opposite cluster,
    %I do not need the calculates of the distances for the same vector. This
    %section deletes those sections not needed.
    for ii = 1:k
        Silhouette(ii,ii).Cluster = [];
    end
    
    %This averages the distance measure calculated for Silhouette.Cluster
    %variable when cluster (k) equals 2. 
    if k == 2
        Silhouette_Between = struct([]);
        for ii = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
                Silhouette_Between(ii).Ave(a,1) = mean(Silhouette(Cluster_k(ii,1),ii).Cluster(a,:));
            end
        end
    end
    
    if k > 2
        Silhouette_Between_Matrix = zeros(k,k);
        for ii = 1:k
            for a = 1:k
                Silhouette_Between_Matrix(a,ii) = a;
            end
        end
        
        for ii = 1:k
            Silhouette_Between_Matrix(ii,ii) = 0;
        end
            
        for ii = 1:k
            for b = 1:k
                if Silhouette_Between_Matrix(b,ii) ~=0
                    for a = 1:size(Cluster_Vectors{1,ii},2)
                    Silhouette_Between(b,ii).Ave(a,1) = mean(Silhouette(b,ii).Cluster(a,:));
                    end
                end
            end
        end
        
        
        SC_Between = struct([]);
        for ii = 1:k
            a = 1;
            for b = 1:k
                if Silhouette_Between_Matrix(b,ii) ~=0     
                    SC_Between(1,ii).Total(:,a) = Silhouette_Between(b,ii).Ave;
                a = a+1;
                end
            end
        end
        
        for ii = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
                SC_Between(1,ii).Ave(a,1) = mean(SC_Between(1,ii).Total(a,:));
            end
        end
        
        for ii = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
                SC_Between(1,ii).Min(a,1) = min(SC_Between(1,ii).Total(a,:));
            end
        end
    end
    
    if k > 2
        for ii = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
            Cluster_Distance_Within(ii).SC(a,1) = (SC_Between(ii).Ave(a,1) - Cluster_Distance_Within(ii).Silhouette_Ave(a,1))/max(SC_Between(ii).Ave(a,1),Cluster_Distance_Within(ii).Silhouette_Ave(a,1));
            end
        end
        
        
        for ii = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
            Cluster_Distance_Within(ii).SC(a,1) = (SC_Between(ii).Min(a,1) - Cluster_Distance_Within(ii).Silhouette_Ave(a,1))/max(SC_Between(ii).Min(a,1),Cluster_Distance_Within(ii).Silhouette_Ave(a,1));
            end
        end
    end
    
    if k == 2
        for ii = 1:k
            for a = 1:size(Cluster_Vectors{1,ii},2)
                Cluster_Distance_Within(ii).SC(a,1) = (Silhouette_Between(ii).Ave(a,1) - Cluster_Distance_Within(ii).Silhouette_Ave(a,1))/max(Silhouette_Between(ii).Ave(a,1),Cluster_Distance_Within(ii).Silhouette_Ave(a,1));
            end
        end
    end
    
    
    %This places all silhouette coefficients in one vector to be averaged
    Silhouette_Total = zeros(size(Data,2),1);
    b = 1;
    for ii = 1:k
        for a = 1:length(Cluster_Distance_Within(ii).SC)
            Silhouette_Total(b) = Cluster_Distance_Within(ii).SC(a,1);
            b = b+1;
        end
    end
    
    Silhouette_Coefficient_Cluster = cell(1,k);
    for ii = 1:k
        Silhouette_Coefficient_Cluster{:,ii} = Cluster_Distance_Within(ii).SC;
    end
    
    %Each time-series silhouette coefficient is averaged to create an overall
    %silhouette score.
    Silhouette_Score = mean(Silhouette_Total);
    
    %Calculates silhouette score for each cluster. I don't need this but I
    %should keep it just in case
    %for ii = 1:k
        %Silhouette_Score_Cluster(1,ii) = mean(Cluster_Distance_Within(ii).SC);
    %end
else
    Silhouette_Score = 0;    
end

end