function null_time_series = Null_TimeSeries_PCA(data, numComponents)
    % NULLTIMESERIESPCA Generates a null time-series using PCA.
    % 
    % Inputs:
    %   data - (NxM matrix) time-series data (N time points, M sensors/axes)
    %   numComponents - (scalar) number of principal components to retain
    %
    % Output:
    %   null_time_series - (NxM matrix) reconstructed null time-series
    
    % Ensure data is zero-mean
    data_mean = mean(data, 1);
    data_centered = data - data_mean;
    
    % Perform PCA using Singular Value Decomposition (SVD)
    [U, S, V] = svd(data_centered, 'econ');
    
    % Retain only the specified number of principal components
    S_null = S;
    S_null(numComponents+1:end, numComponents+1:end) = 0;
    
    % Reconstruct the null time-series
    null_time_series = U * S_null * V' + data_mean;
end