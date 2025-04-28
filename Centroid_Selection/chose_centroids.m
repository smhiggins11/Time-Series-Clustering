function T = chose_centroids(T)

columnNames = T.Properties.VariableNames;
columnNames = setdiff(columnNames, columnNames(end), 'stable');

for zz = 1:size(T,1)
    b = 1;
    if contains("Dunn_Index",columnNames)

        idx = find(strcmp(columnNames,'Dunn_Index'));
        DI = T.(columnNames{idx}){zz}; 
        
        list(:,idx) = DI;

        [DI_numbers, ~, ~] = unique(DI,'first');

        DI_values = sort(DI_numbers,1,"ascend");
        
        ranked_values = DI;
        for i = 1:length(DI_values)
            ranked_values(ranked_values == DI_values(i)) = i;
        end

        ranked(:,b) = ranked_values;

        evals_used{b} = 'Dunn_Index';

        b = b+1;
    end

    if contains("Silhouette_Coefficient",columnNames)
        
        idx = find(strcmp(columnNames,'Silhouette_Coefficient'));
        SC = T.(columnNames{idx}){zz};  
        
        list(:,idx) = SC;

        [SC_numbers, ~, ~] = unique(SC,'first');

        SC_values = sort(SC_numbers,1,"ascend");
        
        ranked_values = SC;
        for i = 1:length(SC_values)
            ranked_values(ranked_values == SC_values(i)) = i;
        end

        ranked(:,b) = ranked_values;
        
        evals_used{b} = 'Silhouette_Coefficient';
        
        b = b+1;
    end


    if contains("Elbow_Method",columnNames)
        idx = find(strcmp(columnNames,'Elbow_Method'));
        WTSS = T.(columnNames{idx}){zz}; 
        
        list(:,idx) = WTSS;

        [WTSS_numbers, ~, ~] = unique(WTSS,'first');

        WTSS_values = sort(WTSS_numbers,1,"descend");
        
        ranked_values = WTSS;
        for i = 1:length(WTSS_values)
            ranked_values(ranked_values == WTSS_values(i)) = i;
        end

        ranked(:,b) = ranked_values;
        
        evals_used{b} = 'Elbow_Method';

        b = b+1;
    end
    
    
    if contains("Gap_Statistic",columnNames) 
        
        idx = find(strcmp(columnNames,'Gap_Statistic'));
        Gap = T.(columnNames{idx}){zz}; 
        
        list(:,idx) = Gap;

        [Gap_numbers, ~, ~] = unique(Gap,'first');

        Gap_values = sort(Gap_numbers,1,"ascend");
        
        ranked_values = Gap;
        for i = 1:length(Gap_values)
            ranked_values(ranked_values == Gap_values(i)) = i;
        end

        ranked(:,b) = ranked_values;
        
        evals_used{b} = 'Gap_Statistic';

        b = b+1;
    end

    if contains("CHI",columnNames) 
        
        idx = find(strcmp(columnNames,'CHI'));
        CHI = T.(columnNames{idx}){zz}; 
        
        list(:,idx) = CHI;

        [CHI_numbers, ~, ~] = unique(CHI,'first');

        CHI_values = sort(CHI_numbers,1,"ascend");
        
        ranked_values = CHI;
        for i = 1:length(CHI_values)
            ranked_values(ranked_values == CHI_values(i)) = i;
        end

        ranked(:,b) = ranked_values;
        
        evals_used{b} = 'CHI';

        b = b+1;
    end

    if contains("DBI",columnNames) 
        
        idx = find(strcmp(columnNames,'DBI'));
        DBI = T.(columnNames{idx}){zz}; 
        
        list(:,idx) = DBI;

        [DBI_numbers, ~, ~] = unique(DBI,'first');

        DBI_values = sort(DBI_numbers,1,"descend");
        
        ranked_values = DBI;
        for i = 1:length(DBI_values)
            ranked_values(ranked_values == DBI_values(i)) = i;
        end

        ranked(:,b) = ranked_values;
        
        evals_used{b} = 'DBI';

        b = b+1;
    end
  
    add_rank = sum(ranked,2);
       
    [~, idx] = max(add_rank);
    
    
    pick_idx(zz,1) = idx;
    picked_row(zz,:) = list(idx,:);
    centroids{zz,1} = T.Centroids{zz}{idx};
end

T = array2table(picked_row,'VariableNames',columnNames);
T.Centroids = centroids;