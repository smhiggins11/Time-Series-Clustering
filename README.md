# Time-Series-Clustering

Overview

This repository contains MATLAB scripts for clustering time-series data and evaluating clustering results. The codebase includes implementations of k-means and hierarchical clustering, along with multiple evaluation metrics such as the elbow method, gap statistic, Dunn index, and silhouette coefficient. These tools aim to objectively determine the optimal number of clusters and the most suitable clustering algorithm for time-series datasets.

Features

Time-Series Clustering: Implements k-means and hierarchical clustering for time-series data.

Clustering Evaluation: Includes multiple evaluation metrics (elbow method, gap statistic, Dunn index, silhouette coefficient) to assess clustering performance.

Automated Selection of Cluster Count: Uses various metrics to determine the best number of clusters.

Visualization Tools: Plots time-series data and clustering results for better interpretability.

Installation

Clone the repository:

    git clone https://github.com/smhiggins11/Time-Series-Clustering

Open MATLAB and navigate to the cloned repository folder.

Usage

Load your time-series dataset into MATLAB.

Time-series data must be equal length

Time-series data must be: (length of time-series x number of time-series)

Run Kmeans.m or HCA_Ward.m (Clustering_Codes folder) to apply clustering.

    [Centroids,index,Vectors,Cluster_n,Initial_Centroids,empty_cluster,not_converged] = Kmeans(Distance_Metric,k,Data,Initial_Centroids)

Inputs for K-means:

Distance_Metric   - Distance metric ('DTW' or 'Euclidean')

k                 - Number of clusters

Data              - Data matrix (time-series data as columns)

Initial_Centroids - (Optional) Predefined initial centroids. If initial centroids is blank, code will automaticly generate random centroids from the dataset

Outputs for K-means:

Centroids         - Final cluster centroids

index             - Indices of data points assigned to each cluster

Vectors           - Time-series for each cluster

Cluster_n         - Number of elements in each cluster

Initial_Centroids - Initial centroid indices

empty_cluster     - Number of times clusters were reinitialized due to empty clusters

not_converged     - Number of times clusters failed to converge within 50 iteration

    [Centroids,index,Cluster_Vectors,Cluster_n,Dendrogram_Total] = HCA_Ward(Distance_Metric,Cluster_range,Data)

HCA_Ward performs Hierarchical Clustering using the Ward method.

Inputs for HCA:

Distance_Metric: Metric to calculate distances ('DTW' or 'Euclidean')

Cluster_range: Desired range of clusters (stopping criteria)

Data: Matrix of time series data (columns are subjects)

Outputs for HCA:

Centroids: The centroids (mean) for each cluster

index: Index of subjects assigned to each cluster

Cluster_Vectors: Time-series for each cluster

Cluster_n: The number of elements in each cluster

Dendrogram_Total: The full hierarchical clustering results

Dependencies

MATLAB (R2021a or later recommended)

Need Signal Processing Toolbox for dtw function

Author

Seth Higgins

For any questions or suggestions, please contact me via GitHub or email.
