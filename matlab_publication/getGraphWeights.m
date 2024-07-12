function [result] = getGraphWeights(traitTable, gene)  
    % Expression data from this database is available upon request
    dbfile = 'GeneExpressions.db';
    gene = {'LRP6','SCD','FASN','ELOVL6','IRS1','FOXO1'};
    next = size(traitTable,2)+1;

    for i=1:size(gene,2)
        connection = sqlite(dbfile);
        query = ['select SampleName, dCT from Expression_',gene{i}];
        dCT = fetch(connection,query);
        close(connection);
        ID = double(cell2mat(dCT(:,1)));
        val = double(cell2mat(dCT(:,2)));
        tmp = [array2table(EpiNr) array2table(val)];

        traitTable(:,next) = array2table(NaN);
        for j=1:size(tmp,1)
            for k=1:size(traitTable,1)
                if(tmp.ID(j) == traitTable.ID(k))
                    traitTable(k,next) = array2table(tmp.val(j));
                    break;
                end
            end
        end
        traitTable.Properties.VariableNames{next} = gene{i};
        next = next +1;

    end

    start = size(traitTable,2) - size(gene,2)+1;
    last = size(traitTable,2);
    result = traitTable(:,[1 start:last]);
end
