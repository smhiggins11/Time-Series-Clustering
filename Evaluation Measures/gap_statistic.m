function GS = gap_statistic(alg,Data,k,Within_Cluster_Total_Sum,Initial_Centroids,Distance_Metric)
%% Calculate Gap Statistic

if contains(alg,"K-means")
    Elog_Wk = zeros(1,100);
    for gap_loop = 1:100
        disp('k = '+string(k)+' N = '+gap_loop)
        surrogate_data = Null_TimeSeries_PCA(Data,2); %create surrogate data gap statistic
        [Gap_Centroid,~,Gap_Vectors,~,~,~,~] = Kmeans(Distance_Metric,k,surrogate_data,Initial_Centroids);
        Gap_WCSS = TWSS(k,Gap_Centroid,Gap_Vectors,Distance_Metric);
        Elog_Wk(1,gap_loop) = log(Gap_WCSS);
        clc
    end
    Elog_Wk_Ave = mean(Elog_Wk);    
    GS = Elog_Wk_Ave - log(Within_Cluster_Total_Sum);
end

if contains(alg,"HCA")
    Elog_Wk = cell(1,k);
    for i = 1:k
        Elog_Wk{i} = zeros(1,100);
    end
    for gap_loop = 1:100
        disp(append('N = ',string(gap_loop)))
        surrogate_data = Null_TimeSeries_PCA(Data,2); %create surrogate data gap statistic
        [Gap_Centroid,~,Gap_Vectors,~,~] = HCA_Ward(Distance_Metric,k,surrogate_data);
        for i = 1:k
            Gap_WCSS = TWSS(i,Gap_Centroid{i,1},Gap_Vectors{i,1},Distance_Metric);
            Elog_Wk{i}(1,gap_loop) = log(Gap_WCSS);
        end
        clc
    end
    for i = 1:k
        Elog_Wk_Ave{i} = mean(Elog_Wk{1});  
    end

    for i = 1:k
        GS(i) = Elog_Wk_Ave{i} - log(Within_Cluster_Total_Sum(i));
    end
end