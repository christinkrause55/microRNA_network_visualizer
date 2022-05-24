function [targets] = callSeedPrediction(IDs,mRNAs, sequence_miR, sequence_mRNA)

    %% Create a miRSeq table based on input IDs
    indx = [];
    for i=1:size(IDs,1)
        f = find(strcmp(string(IDs{i}),table2array(sequence_miR(:,2))));
        indx = [indx f];
    end
    input_miRSeq = sequence_miR(indx.',:);
    
    %% Create a mRSeq table based on input target genes
    indx = [];
    for i=1:size(mRNAs,2)
        f = find(strcmp(string(mRNAs{i}),table2array(sequence_mRNA(:,2))))
        if ~isempty(f)
            indx = [indx f];
        else
            disp(string(mRNAs{i}) + " not found")
        end
    end
    input_mRSeq = sequence_mRNA(indx.',:);

    %% Create empty adjacency matrix
    targets = zeros(size(input_miRSeq,1),size(input_mRSeq,1));
    
    %% Loop
    for i=1:size(input_miRSeq,1)
        for j=1:size(input_mRSeq,1)
            [table_out] = getSeedPrediction(table2array(input_miRSeq(i,1:2)), ...
                convertStringsToChars(table2array(input_mRSeq(j,1))),...
                '',0);
            if(sum(table2array(table_out(:,2))) > 0)
                targets(i,j) = 8;
            elseif(sum(table2array(table_out(:,3))) > 0 || sum(table2array(table_out(:,4))) > 0 || sum(table2array(table_out(:,5))) > 0 || sum(table2array(table_out(:,6))) > 0 )
                targets(i,j) = 7;
            end
        end
    end
end

