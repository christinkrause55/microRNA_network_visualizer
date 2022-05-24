function [table_out] = getSeedPrediction(miRSeq, mRNASeq,nameTable,writeTo)
    %miRSeq = ["UGUCAGUUUGUCAAAUACCCCA" "mmu-miR-223-3p"];
    mRNASeq = erase(mRNASeq,' ');
    
%% Init table    
    %Name_miR | 8mer | 7mer_m8 | 7mer_A1 | 6mer | offset_6mer | %
    %AU_content | in_proximity_same (<34 nt) | in_proximity_other (<8-40 nt)|%
    %3'-miR-Binding at 12 -17 nt | > 15 nt Stop codon | relative location %
    table_out = table(miRSeq(1,2),0,0,0,0,0,0,false,convertCharsToStrings(' '),convertCharsToStrings(' '),false,0.0);
    table_out.Properties.VariableNames{1} = 'miRNA';
    table_out.Properties.VariableNames{2} = 'Pos_8mer';
    table_out.Properties.VariableNames{3} = 'Pos_7mer_m8';
    table_out.Properties.VariableNames{4} = 'Pos_7mer_A1';
    table_out.Properties.VariableNames{5} = 'Pos_6mer';
    table_out.Properties.VariableNames{6} = 'Pos_offset_6mer';
    table_out.Properties.VariableNames{7} = 'AU_content_30nt';
    table_out.Properties.VariableNames{8} = 'In_proximity_same';
    table_out.Properties.VariableNames{9} = 'In_proximity_other';
    table_out.Properties.VariableNames{10} = '3_miR_binding';
    table_out.Properties.VariableNames{11} = 'Stop_codon';
    table_out.Properties.VariableNames{12} = 'relative_within_UTR';
    table_idx = 1;
    
%% Loop for miRNA sequences -> process_seed, find seeds in mRNA
    for seq=1:size(miRSeq)
        disp(seq)
        seed_miR = process_seed(miRSeq(seq,1));
        
        % Fill array for 8mer
        positions = contains_8mer(seed_miR, mRNASeq);
        for j=1:size(positions,2)
            table_out(table_idx,1) = array2table(miRSeq(seq,2));
            table_out(table_idx,2) = array2table(positions(j));
            table_out(table_idx,7) = array2table(containsAU(positions(j), mRNASeq));
            [a,b] = isBindingAlsoUpstream(positions(j),miRSeq(seq,1),mRNASeq);
            table_out(table_idx,10) = array2table(convertCharsToStrings(b));
            table_out(table_idx,11) = array2table(isBehindStop(positions(j),true));
            table_out(table_idx,12) = array2table(whereInUTR(positions(j),mRNASeq));
            
            table_idx = table_idx +1;
        end
        
        % Fill array for 7mer m8
        positions = contains_7mer_m8(seed_miR, mRNASeq);
        for j=1:size(positions,2)
            table_out(table_idx,1) = array2table(miRSeq(seq,2));
            table_out(table_idx,3) = array2table(positions(j));
            table_out(table_idx,7) = array2table(containsAU(positions(j), mRNASeq));
            [a,b] = isBindingAlsoUpstream(positions(j),miRSeq(seq,1),mRNASeq);
            table_out(table_idx,10) = array2table(convertCharsToStrings(b));
            table_out(table_idx,11) = array2table(isBehindStop(positions(j),true));
            table_out(table_idx,12) = array2table(whereInUTR(positions(j),mRNASeq));
            
            table_idx = table_idx +1;
        end
        
        % Fill array for 7mer A1
        positions = contains_7mer_A1(seed_miR, mRNASeq);
        for j=1:size(positions,2)
            table_out(table_idx,1) = array2table(miRSeq(seq,2));
            table_out(table_idx,4) = array2table(positions(j));
            table_out(table_idx,7) = array2table(containsAU(positions(j), mRNASeq));
            [a,b] = isBindingAlsoUpstream(positions(j),miRSeq(seq,1),mRNASeq);
            table_out(table_idx,10) = array2table(convertCharsToStrings(b));
            table_out(table_idx,11) = array2table(isBehindStop(positions(j),true));
            table_out(table_idx,12) = array2table(whereInUTR(positions(j),mRNASeq));
            
            table_idx = table_idx +1;
        end
        
        positions = contains_6mer(seed_miR, mRNASeq);
        for j=1:size(positions,2)
            table_out(table_idx,1) = array2table(miRSeq(seq,2));
            table_out(table_idx,5) = array2table(positions(j));
            table_out(table_idx,7) = array2table(containsAU(positions(j), mRNASeq));
            [a,b] = isBindingAlsoUpstream(positions(j),miRSeq(seq,1),mRNASeq);
            table_out(table_idx,10) = array2table(convertCharsToStrings(b));
            table_out(table_idx,11) = array2table(isBehindStop(positions(j),true));
            table_out(table_idx,12) = array2table(whereInUTR(positions(j),mRNASeq));
            
            table_idx = table_idx +1;
        end
        
        positions = contains_offset_6mer(seed_miR, mRNASeq);
        for j=1:size(positions,2)
            table_out(table_idx,1) = array2table(miRSeq(seq,2));
            table_out(table_idx,6) = array2table(positions(j));
            table_out(table_idx,7) = array2table(containsAU(positions(j), mRNASeq));
            [a,b] = isBindingAlsoUpstream(positions(j),miRSeq(seq,1),mRNASeq);
            table_out(table_idx,10) = array2table(convertCharsToStrings(b));
            table_out(table_idx,11) = array2table(isBehindStop(positions(j),true));
            table_out(table_idx,12) = array2table(whereInUTR(positions(j),mRNASeq));
            
            table_idx = table_idx +1;
        end
        
        % Compare for same miRNA within region of 34 nt?
        %miR_idx = find(strcmp(table2array(table_out(:,1)),convertCharsToStrings(miRSeq(seq,2))));
        
    end
    
    % in proximity others
    for m=1:size(table_out,1)
        position_idx = zeros(1,size(table_out,2));
        out_append_names = "";
        for n =2:6
            position_idx(1,n) = table2array(table_out(m,n)) ~= 0; 
        end
        value_pos = table2array(table_out(m,find(position_idx==1)));
        in_prox = ((table2array(table_out(:,2)) <= (value_pos+34) & table2array(table_out(:,2)) >= (value_pos-34)) | (table2array(table_out(:,3)) <= (value_pos+34) & table2array(table_out(:,3)) >= (value_pos-34)) | (table2array(table_out(:,4)) <= (value_pos+34) & table2array(table_out(:,4)) >= (value_pos-34)) | (table2array(table_out(:,5)) <= (value_pos+34) & table2array(table_out(:,5)) >= (value_pos-34)) | (table2array(table_out(:,6)) <= (value_pos+34) & table2array(table_out(:,6)) >= (value_pos-34))) & ~strcmp(table2array(table_out(:,1)),table2array(table_out(m,1)));
        
        if(sum(in_prox) ~= 0)
            names_miR = table_out(in_prox(:,1),1);
            out_names = unique(names_miR);
            for o=1:size(names_miR,1)
                out_append_names = append(out_append_names,table2array(names_miR(o,1)),'; ');
            end    
        else
            out_append_names = "-";
        end
        table_out(m,9) = array2table(out_append_names);
    end
    
    
%% Ploting
if(false)
    subplot(3,2,1)
    values = table2array(table_out((table2array(table_out(:,2))~= 0),2));
    histogram(values,round(size(mRNASeq,2)))
    title('8mer')
    ylabel('Count')
    xlabel('nt of target gene')
    subplot(3,2,2)
    values = table2array(table_out((table2array(table_out(:,3))~= 0),3));
    histogram(values,round(size(mRNASeq,2)))
    title('7mer m8')
    ylabel('Count')
    xlabel('nt of target gene')
    subplot(3,2,3)
    values = table2array(table_out((table2array(table_out(:,4))~= 0),4));
    histogram(values,round(size(mRNASeq,2)))
    title('7mer A1')
    ylabel('Count')
    xlabel('nt of target gene')
    subplot(3,2,4)
    values = table2array(table_out((table2array(table_out(:,5))~= 0),5));
    histogram(values,round(size(mRNASeq,2)))
    title('6mer')
    ylabel('Count')
    xlabel('nt of target gene')
    subplot(3,2,5)
    values = table2array(table_out((table2array(table_out(:,6))~= 0),6));
    histogram(values,round(size(mRNASeq,2)))
    title('offset 6mer')
    ylabel('Count')
    xlabel('nt of target gene')
    subplot(3,2,6)
    values = table2array(table_out((table2array(table_out(:,2))~= 0 & ~strcmp(table2array(table_out(:,10)),"0 - 0")),2));
    histogram(values,round(size(mRNASeq,2)))
    title('8mer and additional seed')
    ylabel('Count')
    xlabel('nt of target gene')
end
%% Additional functions for loops
   
% Get 7nt seed of miRNA on 5' end
    function seed = process_seed(miRSeq)
        seed_prep = extractBetween(miRSeq,2,8);
        % From RNA into DNA
        seed_prep = strrep(seed_prep,'U','T');
        % Build reverse and complement for sequence to find on mRNA
        seed_prep = reverse(seed_prep);
        temp = char(seed_prep(1));
        seed = ' ';
        for s = 1:size(temp,2)
           if temp(s) == 'A'
               seed = strcat(seed,'T');
           elseif temp(s) == 'T'
               seed = strcat(seed,'A');
           elseif temp(s) == 'G'
               seed = strcat(seed,'C');
           else 
               seed = strcat(seed,'G');
           end
        end
    end

% Get 6 nt seed (12-17) of miRNA on 3' end
    function seed = upstream_seed(miRSeq)
        %Extraction, conversion and reverse
        if(size(char(miRSeq),2)>= 17)
            seed_prep = reverse(strrep(extractBetween(miRSeq,12,17),'U','T'));
        else
            seed_prep = reverse(strrep(extractBetween(miRSeq,12,size(char(miRSeq),2)),'U','T'));
        end
        temp = char(seed_prep(1));
        seed = ' ';
        for s = 1:size(temp,2)
           if temp(s) == 'A'
               seed = strcat(seed,'T');
           elseif temp(s) == 'T'
               seed = strcat(seed,'A');
           elseif temp(s) == 'G'
               seed = strcat(seed,'C');
           else 
               seed = strcat(seed,'G');
           end
        end
    end

% Get all positions of 8mer (7 nt complementary binding + A opposite of nt 1)
    function [k] = contains_8mer (seed, mRNASeq)
        k = strfind(mRNASeq, strcat(seed,'A'));
    end

% Get all positions of 7mer_A1 (6 nt complementary biding + A opposite of nt 1)
    function [k] = contains_7mer_A1 (seed, mRNASeq)
        k = strfind(mRNASeq, strcat(seed(2:end),'A'));
    end

% Get all positions of 7mer_m8 (7 nt complementary binding)
    function [k] = contains_7mer_m8 (seed, mRNASeq)
        k = strfind(mRNASeq, seed);
    end

% Get all positions of 6mer (6 nt complementary binding, nt 2-7)
    function [k] = contains_6mer (seed, mRNASeq)
        k = strfind(mRNASeq, strcat(seed(2:end)));
    end

% Get all positions of offset 6mer (6 nt complementary binding, nt 3-8)
    function [k] = contains_offset_6mer (seed, mRNASeq)
        k = strfind(mRNASeq, strcat(seed(1:end-1)));
    end

% Get boolean whether opstream of 5' seed there is an additional 3' seed and its
% respective positions
    function [contains,temp] = isBindingAlsoUpstream(pos,miRSeq,mRNASeq)
        contains = false;
        seed = upstream_seed(miRSeq);
        limit_upstream_seed = pos - 5 - 6;
        position_secondSeed = zeros(1,2);
        idx = 1;
        
        if(limit_upstream_seed >= 1)
            roi = extractBetween(mRNASeq,limit_upstream_seed,pos-1);
        else
            roi = extractBetween(mRNASeq,1,pos-1);
        end
        
        for i=1:2
            k = strfind(roi,seed(i:end));
           
            if(~isempty(k{1}))
                contains = true;
                position_secondSeed(idx,1) = k{1}(1) + limit_upstream_seed-1;
                position_secondSeed(idx,2) = size(seed(i:end),2);
                idx = idx+1;
            end
            k = strfind(roi,seed(i:end-1));
            if(~isempty(k{1}))
                contains = true;
                position_secondSeed(idx,1) = k{1}(1) + limit_upstream_seed-1;
                position_secondSeed(idx,2) = size(seed(i:end),2);
                idx = idx+1;
            end
        end
        % Substract 1 from upper limit, because new 1 of substring is upper
        % limit + relative position within new substring
        position_out = unique(position_secondSeed,'rows');
        
        temp = append(num2str(position_out(1,1)),' - ', num2str(position_out(1,2)));
        if(size(position_out) > 1)
            for p=2:size(position_out,2)
                temp = append(temp, ';',num2str(position_out(p,1)),' - ', num2str(position_out(p,2)));
            end
        end
        
            
    end

% Counts the percentage of AU 30 nt around seed
    function [percent_AU] = containsAU(pos, mRNASeq)
        limit_seed_dw = pos - 30;
        limit_seed_up = pos + 7 + 30;
        if(limit_seed_up > size(mRNASeq,2))
            limit_seed_up = size(mRNASeq,2);
        end
        if(limit_seed_dw >= 1 && limit_seed_up <= size(mRNASeq,2))
            roi = extractBetween(mRNASeq,limit_seed_dw,limit_seed_up);
        elseif(limit_seed_dw <= 0 && limit_seed_up <= size(mRNASeq,2))
            roi = extractBetween(mRNASeq,1,limit_seed_up);
        elseif(limit_seed_dw >= 1 && limit_seed_up > size(mRNASeq,2))
            roi = extractBetween(mRNASeq,limit_seed_dw,size(mRNASeq,2)+limit_seed_up);
        else
            roi = extractBetween(mRNASeq,1,size(mRNASeq,2));
        end

        counter_A_T = 0;
        roi_iteratorString = roi{1};

        for i=1:size(roi_iteratorString,2)
            if(roi_iteratorString(i) == 'A' || roi_iteratorString(i) == 'T')
                counter_A_T = counter_A_T +1;
            end
        end
        percent_AU = counter_A_T/size(roi_iteratorString,2);
        
    end

% Needs a flag, whether the mRNA sequence was only 3'UTR
    function [behindStop] = isBehindStop(pos,flag)
        behindStop = false;
        if(flag)
            behindStop = pos > 15;
        end
    end 

% Calculates the relative location of a position within an UTR and
% calculates a simple weighting score. Zero for seed within first 15 nt of
% UTR.
    function [relative_inUTR, isEffective] = whereInUTR(pos,mRNASeq)
        sizeUTR = size(mRNASeq,2);
        relative_inUTR = pos/sizeUTR;
        noEffect = 15/sizeUTR;
        
        if(relative_inUTR > noEffect)
            if(relative_inUTR <= 0.25 || relative_inUTR >= 0.75)
                isEffective = 2;
            elseif (relative_inUTR >= 0.6 || relative_inUTR <= 0.4)
                isEffective = 1.5;
            else 
                isEffective = 1;
            end
        else
            isEffective = 0;
        end
    end
    
if(writeTo == 1)
    writetable(table_out,nameTable,'WriteVariableNames',true)
end

end

