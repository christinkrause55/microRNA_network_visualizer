function [result] = getGraphWeights_mouseLRGs_VSG(gene, animals)    
    dbfile = 'GeneExpressions.db';
    %next = size(traitTable,2)+1;

    connection = sqlite(dbfile);
    query = ['select * from VSG_LRG_expression'];
    dCT = fetch(connection,query);
    close(connection);
    
    result = zeros(sum(animals),size(gene,2));
    
    for i=1:size(gene,2)
        indx = find(strcmp(gene{i},dCT(:,1)));
        tmp = cell2mat(dCT(indx,[2:sum(animals)+1]));
        result(:,i) = tmp;
    end

    result = array2table(result);
    for i=1:size(gene,2)
        result.Properties.VariableNames{i} = gene{i};
    end
end