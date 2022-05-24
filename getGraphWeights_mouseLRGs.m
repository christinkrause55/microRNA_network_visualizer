function [result] = getGraphWeights_mouseLRGs(gene, animals, protein)    
    dbfile = 'GeneExpressions.db';
    %next = size(traitTable,2)+1;

    connection = sqlite(dbfile);
    if(~protein)
        query = ['select * from DZD_LRG_expression'];
        dCT = fetch(connection,query);
        close(connection);
    else
        query = ['select * from DZD_proteom'];
        dCT = fetch(connection,query);
        close(connection);
    end
    
    result = zeros(sum(animals),size(gene,2));
    
    for i=1:size(gene,2)
        indx = find(strcmp(gene{i},dCT(:,1)));
        tmp = cell2mat(dCT(indx,[2:sum(animals)+1]));
        result(:,i) = tmp;
        indx2 = result(:,i) == 0;
        result(indx2,i) = NaN;
    end

    result = array2table(result);
    for i=1:size(gene,2)
        result.Properties.VariableNames{i} = gene{i};
    end
end