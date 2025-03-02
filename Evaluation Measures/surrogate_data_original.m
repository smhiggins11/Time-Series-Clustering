function surrogate_data = surrogate_data_original(Matrix)

Matrix_Sum = zeros(length(Matrix),1);

for ii = 1:size(Matrix,2)
    Matrix_Sum = Matrix_Sum+Matrix(:,1);
end

Matrix_Ave = Matrix_Sum/size(Matrix,2);

for ii = 1:size(Matrix,2)
    Matrix_Deviation(:,ii) = (Matrix(:,ii)-Matrix_Ave).^2;
end

Matrix_Deviation_Sum = zeros(length(Matrix),1);
for ii = 1:size(Matrix,2)
    Matrix_Deviation_Sum = Matrix_Deviation_Sum+Matrix_Deviation(:,ii);
end

Matrix_Deviation_Ave = Matrix_Deviation_Sum/size(Matrix,2);

Matrix_SD_plus = sqrt(Matrix_Deviation_Ave);
Matrix_SD_plus_Ave = mean(Matrix_SD_plus);

Matrix_SD_min = Matrix_SD_plus*-1;
Matrix_SD_min_Ave = mean(Matrix_SD_min);

for ii = 1:size(Matrix,2)
    surrogate_data(:,ii) = Matrix(:,ii)+(Matrix_SD_plus_Ave + (Matrix_SD_min_Ave-Matrix_SD_plus_Ave) .* rand(1,1));
end
