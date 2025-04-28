function centroids = initializeCentroids(Method,Data, k)
    % Inputs:
    % - timeSeriesData: Cell array of time-series (each cell contains a 1D array)
    % - k: Number of clusters
    % Output:
    % - centroids: k initial centroids selected using K-means++ with DTW
    
    lengthseries = size(Data,1);
    numSeries = size(Data,2);
    centroids = zeros(lengthseries,k);

    % Step 1: Select the first centroid 
    % 
    % There are six methods I could use
    
    % First Method is to choose all centroids randomly
    if contains(Method, 'Random-All')
        Initial_Centroids = randperm(size(Data,2), k);
        centroids = Data(:, Initial_Centroids);
        return
    end

    % Second Method is to choose first centroid randomly
    if contains(Method,"KMeans++")
        firstIdx = randi(numSeries);
        centroids(:,1) = Data(:,firstIdx);
    end

    % Third method is to choose first centroid as time-series nearest to
    % the global mean
    
    %this averages time-series taking into account if time-series are different
    %lengths.
    if contains(Method,"MeanDTW++")
        TS_mean = zeros(size(Data,1),1); %preallocate variable
        for i = 1:size(Data,1)
            len = sum(~isnan(Data(i,:)),'all');
            summation = nansum(Data(i,:)); 
            TS_mean(i,1) = summation/len;
        end
    
        for i = 1:numSeries
            dist(i) = dtw(TS_mean(~isnan(TS_mean)),Data(~isnan(Data(:,i)),i));
        end
        
        [min_dist,min_idx] = min(dist);
        centroids(:,1) = Data(:,min_idx);
    end

    % Fourth method is to choose first centroid as time-series that is
    % furthest away from the mean
    if contains(Method,"MaxVar++")
        TS_mean = zeros(size(Data,1),1); %preallocate variable
        for i = 1:size(Data,1)
            len = sum(~isnan(Data(i,:)),'all');
            summation = nansum(Data(i,:)); 
            TS_mean(i,1) = summation/len;
        end
    
        for i = 1:numSeries
            dist(i) = dtw(TS_mean(~isnan(TS_mean)),Data(~isnan(Data(:,i)),i));
        end
        [max_dist,max_idx] = max(dist);
        centroids(:,1) = Data(:,max_idx);
    end
    
    % Fifth method chooses first centroid using HCA. First centroid will be
    % selected from the largest cluster and will be the closest time-series
    % to the centroid.
    if contains(Method,"HCA++")
        [Centroids,~,Cluster_Vectors,Cluster_n,~] = HCA_Ward("DTW",k,Data);
        largestCluster = mode(Cluster_n{k,1});
        largestCluster_idx = find(Cluster_n{k,1} == largestCluster);
        if length(largestCluster_idx) > 1
            largestCluster_idx = largestCluster_idx(1,1);
        end
        for i = 1:size(Cluster_Vectors{k,1}{1,largestCluster_idx},2)
            dist(i) = dtw( ...
                Centroids{k,1}(~isnan(Centroids{k,1}(:,largestCluster_idx)), ...
                largestCluster_idx),Cluster_Vectors{k,1}{1,largestCluster_idx}(~isnan(Cluster_Vectors{k,1}{1,largestCluster_idx}(:,i)),i));
        end
        min_dist = min(dist);
        centroid_idx = dist == min_dist;
        centroids(:,1) = Data(:,centroid_idx);
    end
    
    % Sixth method chooses all centroids using HCA
    if contains(Method,'Full-HCA')
        [Centroids,~,Cluster_Vectors,Cluster_n,~] = HCA_Ward("DTW",k,Data);
        centroids = cell(k,1);
        for iii = 1:k
            for ii = 1:size(Cluster_n{iii},2)
                distance = zeros(1,Cluster_n{iii}(ii));
                for i = 1:Cluster_n{iii}(ii)
                    distance(i) = dtw(Cluster_Vectors{iii,1}{1,ii}(:,i),Centroids{iii,1}(:,ii));
                end
                indx(:,ii) = Cluster_Vectors{iii,1}{1,ii}(:,distance == min(distance));
            end
            centroids{iii} = indx;
        end
        return
    end

    % Step 2: Select remaining centroids using DTW-based probability
    % This section will be used for each method except for Random-All and
    % Full HCA methods
    for j = 2:k
        dists = zeros(1, numSeries); % Store min DTW distances for each series
        
        % Compute DTW distance to the nearest selected centroid
        for i = 1:numSeries
            minDist = inf;
            for c = 1:j-1
                dist = dtw(Data(~isnan(Data(:,i)),i), centroids(~isnan(centroids(:,c)),c));
                minDist = min(minDist, dist);
            end
            dists(i) = minDist^2; % Use squared distance for probability weighting
        end

        % Select the next centroid with probability proportional to dists^2
        probDist = dists / sum(dists);
        selectedIdx = randsample(1:numSeries, 1, true, probDist);
        centroids(:,j) = Data(:,selectedIdx);
    end
end