function GS = gap_statistic(Distance_Metric,Data,k,Cluster_Centroids,Cluster_Vectors,Initial_Centroids)
%% Calculate Gap Statistic

Within_Cluster_Total_Sum = TWSS(k,Cluster_Centroids,Cluster_Vectors,Distance_Metric);

Elog_Wk = zeros(1,100);
for gap_loop = 1:100
    disp('k = '+string(k)+' N = '+gap_loop)
    surrogate_data = surrogate_data_original(Data); %create surrogate data gap statistic
    [Gap_Centroid,~,Gaps_Vectors,~,~,~,~] = Kmeans(Distance_Metric,k,surrogate_data,Initial_Centroids);
    Gap_WCSS = TWSS(k,Gap_Centroid,Gaps_Vectors,Distance_Metric);
    Elog_Wk(1,gap_loop) = log(Gap_WCSS);
    clc
end
Elog_Wk_Ave = mean(Elog_Wk);    
GS = Elog_Wk_Ave - log(Within_Cluster_Total_Sum);