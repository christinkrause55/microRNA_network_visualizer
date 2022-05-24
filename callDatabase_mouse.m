dbfile = 'GeneExpressions.db';

%% Get miRNA Names
s='''';
termMIR = strcat(s,'mmu-%',s);
connection = sqlite(dbfile);
query = ['select ID,TranscriptID from DZD_miR4_values where DZD_miR4_values.TranscriptID like ', ...
    termMIR];
miRNA_IDs = fetch(connection,query);
close(connection);

table_pvalues = array2table(zeros(size(miRNA_IDs,1),46));
table_pvalues.Properties.VariableNames{1} = 'pHFD_chow';
table_pvalues.Properties.VariableNames{2} = 'eHFD_chow';
table_pvalues.Properties.VariableNames{3} = 'pHC_chow';
table_pvalues.Properties.VariableNames{4} = 'eHC_chow';
table_pvalues.Properties.VariableNames{5} = 'isexpressed';
table_pvalues.Properties.VariableNames{6} = 'meanChow';
table_pvalues.Properties.VariableNames{7} = 'pCR_chow';
table_pvalues.Properties.VariableNames{8} = 'eCR_chow';
table_pvalues.Properties.VariableNames{9} = 'pEx_chow';
table_pvalues.Properties.VariableNames{10} = 'eEx_chow';
table_pvalues.Properties.VariableNames{11} = 'pVSG_d9';
table_pvalues.Properties.VariableNames{12} = 'eVSG_d9';
table_pvalues.Properties.VariableNames{13} = 'pVSG_d35';
table_pvalues.Properties.VariableNames{14} = 'eVSG_d35';
table_pvalues.Properties.VariableNames{15} = 'pSham_d9_d35';
table_pvalues.Properties.VariableNames{16} = 'eSham_d9_d35';
table_pvalues.Properties.VariableNames{17} = 'pVSG_d9_d35';
table_pvalues.Properties.VariableNames{18} = 'eVSG_d9_d35';
table_pvalues.Properties.VariableNames{19} = 'pNAS_human';
table_pvalues.Properties.VariableNames{20} = 'eNAS_human';
table_pvalues.Properties.VariableNames{21} = 'pT2D_human';
table_pvalues.Properties.VariableNames{22} = 'eT2D_human';
table_pvalues.Properties.VariableNames{23} = 'pBMI_human';
table_pvalues.Properties.VariableNames{24} = 'eBMI_human';
table_pvalues.Properties.VariableNames{25} = 'isexpressed_human';
table_pvalues.Properties.VariableNames{26} = 'meanAll_human';
table_pvalues.Properties.VariableNames{27} = 'pHepatoLipid_human';
table_pvalues.Properties.VariableNames{28} = 'eHepatoLipid_human';
table_pvalues.Properties.VariableNames{29} = 'pGlc';
table_pvalues.Properties.VariableNames{30} = 'eGlc';
table_pvalues.Properties.VariableNames{31} = 'pHbA1c';
table_pvalues.Properties.VariableNames{32} = 'eHbA1c';
table_pvalues.Properties.VariableNames{33} = 'pTriglycerides';
table_pvalues.Properties.VariableNames{34} = 'eTriglycerides';
table_pvalues.Properties.VariableNames{35} = 'pInsulin';
table_pvalues.Properties.VariableNames{36} = 'eInsulin';
table_pvalues.Properties.VariableNames{37} = 'pCpep';
table_pvalues.Properties.VariableNames{38} = 'eCpep';
table_pvalues.Properties.VariableNames{39} = 'pHOMAIR';
table_pvalues.Properties.VariableNames{40} = 'eHOMAIR';

table_pvalues.Properties.VariableNames{41} = 'pCR_HFD';
table_pvalues.Properties.VariableNames{42} = 'eCR_HFD';
table_pvalues.Properties.VariableNames{43} = 'pHC_HFD';
table_pvalues.Properties.VariableNames{44} = 'eHC_HFD';
table_pvalues.Properties.VariableNames{45} = 'pEx_HFD';
table_pvalues.Properties.VariableNames{46} = 'eEx_HFD';

table_pvalues(:,47) = array2table(miRNA_IDs(:,2));
table_pvalues.Properties.VariableNames{47} = 'TranscriptID';

for k=1:size(miRNA_IDs,1)
    disp('moving to position ' + string(k))
    
    %%% MOUSE DZD 
    %% Get miRNA Information
    connection = sqlite(dbfile);
    term = strcat(s,miRNA_IDs{k},s);
    term2 = strcat(s,'sample',s);
    query = ['select * from DZD_miR4_values where DZD_miR4_values.ID like ', ...
        term, ' or DZD_miR4_values.TranscriptID = ',...
        term2];

    result = fetch(connection,query);
    close(connection);

    col1 = array2table(cell2mat(result(:,3:end)).');
    col1.Properties.VariableNames{1} = 'MouseID';
    col1.Properties.VariableNames{2} = miRNA_IDs{k,2};
    
    Chow = (1:6).';
    HFD = (7:12).';
    HC = (13:18).';
    CR = (19:24).';
    Ex = (25:30).';
    
    % Row 1:2 Chow vs HFD
    data = [col1(Chow,2); col1(HFD,2)];
    data(Chow,2) = array2table(1);
    mdl = fitlm(data);
    table_pvalues(k,1) = mdl.Coefficients(2,4);
    table_pvalues(k,2) = array2table(mean(table2array(col1(HFD,2))) /mean(table2array(col1(Chow,2))));
    % Row 3:4 HC vs Chow
    data = [col1(Chow,2); col1(HC,2)];
    data(Chow,2) = array2table(1);
    mdl = fitlm(data);
    table_pvalues(k,3) = mdl.Coefficients(2,4);
    table_pvalues(k,4) = array2table(mean(table2array(col1(HC,2))) /mean(table2array(col1(Chow,2))));
    % Row 6 MeanChow
    table_pvalues(k,6) = array2table(mean(table2array(data(Chow,1))));
    % Row 7:8 CR
    data = [col1(Chow,2); col1(CR,2)];
    data(Chow,2) = array2table(1);
    mdl = fitlm(data);
    table_pvalues(k,7) = mdl.Coefficients(2,4);
    table_pvalues(k,8) = array2table(mean(table2array(col1(CR,2))) /mean(table2array(col1(Chow,2))));
    % Row 9:10 Ex
    data = [col1(Chow,2); col1(Ex,2)];
    data(Chow,2) = array2table(1);
    mdl = fitlm(data);
    table_pvalues(k,9) = mdl.Coefficients(2,4);
    table_pvalues(k,10) = array2table(mean(table2array(col1(Ex,2))) /mean(table2array(col1(Chow,2))));    
    % Row 41:42 CR vs HFD
    data = [col1(HFD,2); col1(CR,2)];
    data(HFD,2) = array2table(1);
    mdl = fitlm(data);
    table_pvalues(k,41) = mdl.Coefficients(2,4);
    table_pvalues(k,42) = array2table(mean(table2array(col1(CR,2))) /mean(table2array(col1(HFD,2))));  
    % Row 43:44 HC vs HFD
    data = [col1(HFD,2); col1(HC,2)];
    data(HFD,2) = array2table(1);
    mdl = fitlm(data);
    table_pvalues(k,43) = mdl.Coefficients(2,4);
    table_pvalues(k,44) = array2table(mean(table2array(col1(HC,2))) /mean(table2array(col1(HFD,2))));  
    % Row 41:42 Ex vs HFD
    data = [col1(HFD,2); col1(Ex,2)];
    data(HFD,2) = array2table(1);
    mdl = fitlm(data);
    table_pvalues(k,45) = mdl.Coefficients(2,4);
    table_pvalues(k,46) = array2table(mean(table2array(col1(Ex,2))) /mean(table2array(col1(HFD,2))));  
    
    
    %%% MOUSE VSG 
    %% Get miRNA Information
    connection = sqlite(dbfile);
    term3 = strcat(s,'Group2',s);
    query = ['select * from VSG_miR4_values where VSG_miR4_values.ID like ', ...
        term, ' or VSG_miR4_values.TranscriptID = ',...
        term2, ' or VSG_miR4_values.TranscriptID = ', term3];

    result = fetch(connection,query);
    close(connection);

    col1 = array2table(cell2mat(result(:,3:end)).');
    col1.Properties.VariableNames{1} = 'MouseID';
    col1.Properties.VariableNames{2} = miRNA_IDs{k,2};
    col1.Properties.VariableNames{3} = 'Group';
    
    Sham_d9 = table2array(col1(:,3)) == 1;
    VSG_d9 = table2array(col1(:,3)) == 2;
    Sham_d35 = table2array(col1(:,3)) == 3;
    VSG_d35 = table2array(col1(:,3)) == 4;
    
    % Row 11:12 VSG_d9
    data = [col1(Sham_d9,2:3); col1(VSG_d9,2:3)];
    mdl = fitlm(data);
    table_pvalues(k,11) = mdl.Coefficients(2,4);
    table_pvalues(k,12) = array2table(mean(table2array(col1(VSG_d9,2))) /mean(table2array(col1(Sham_d9,2))));
    % Row 13:14 VSG_d35
    data = [col1(Sham_d35,2:3); col1(VSG_d35,2:3)];
    mdl = fitlm(data);
    table_pvalues(k,13) = mdl.Coefficients(2,4);
    table_pvalues(k,14) = array2table(mean(table2array(col1(VSG_d35,2))) /mean(table2array(col1(Sham_d35,2))));
    % Row 15:16 Sham_d9_d35
    data = [col1(Sham_d35,2:3); col1(Sham_d9,2:3)];
    mdl = fitlm(data);
    table_pvalues(k,15) = mdl.Coefficients(2,4);
    table_pvalues(k,16) = array2table(mean(table2array(col1(Sham_d35,2))) /mean(table2array(col1(Sham_d9,2)))); 
    % Row 17:18 VSG_d9_d35
    data = [col1(VSG_d35,2:3); col1(VSG_d9,2:3)];
    mdl = fitlm(data);
    table_pvalues(k,17) = mdl.Coefficients(2,4);
    table_pvalues(k,18) = array2table(mean(table2array(col1(VSG_d35,2))) /mean(table2array(col1(VSG_d9,2))));  
    
     %% MOUSE: Is it expressed?
     connection = sqlite('miRNA.db');
     query3 = ['select mmuLiver_intervention_miRNA4.ExpressedinAllConditions ',...
         ' from mmuLiver_intervention_miRNA4 ',...
         ' where mmuLiver_intervention_miRNA4.ID = ', term];
     result = fetch(connection,query3);
     close(connection);
     if(result{1}=='F')
         table_pvalues(k,5) = array2table(0);
     else
         table_pvalues(k,5) = array2table(1);
     end
    
    %%% HUMAN
    connection = sqlite(dbfile);
    term = strcat(s,replace(miRNA_IDs{k,2},'mmu','hsa'),s);
    term2 = strcat(s,'sample',s);

    query = ['select * from values_miR4_col where values_miR4_col.TranscriptID like ', ...
        term, ' or values_miR4_col.TranscriptID = ',...
        term2];

    result = fetch(connection,query);
    close(connection);
    if (size(result,1) > 1)
        col1 = array2table(cell2mat(result(:,3:end)).');
        col1.Properties.VariableNames{1} = 'EpiNr';
        col1.Properties.VariableNames{2} = miRNA_IDs{k,2};

        %% Get Parameters
        query2 = ['select UKE_Cohort_Parameters.Epi_Nr, T2DnachUKE, Geschlecht_2_w_1_m, OPBMI, Age, OPNASHScore, HbA1c, ISGlucose, ISTriglyceride, Insulin_pmol_L_525_max, HOMA_IR, Cpeptid_pmol_L, Verfettung_Hepatozyten, OPFibrose ',...
            'from UKE_Cohort_Parameters, UKE_Cohort_Scores, UKE_Cohort_FattyContent ',...
            'where UKE_Cohort_Parameters.Epi_Nr = UKE_Cohort_Scores.Epi_Nr and UKE_Cohort_Parameters.Epi_Nr = UKE_Cohort_FattyContent.Epi_Nr'];
        connection = sqlite(dbfile);
        result = fetch(connection,query2);
        close(connection);
        parameters = {'EpiNr', 'T2DnachUKE', 'Geschlecht_2_w_1_m', 'OPBMI', 'Age', 'OPNASHScore','HbA1c','ISGlucose','ISTriglyceride','Insulin','HOMA_IR', 'Cpeptid', 'Hepatolipid','Fibrose'};

        col2 = array2table(zeros(size(result)));
        for i=1:1
            col2(:,i) = array2table(cell2mat(result(:,i)));
            col2.Properties.VariableNames{i} = parameters{i};
        end

        for i=2:14
            col2(:,i) = array2table(str2double(result(:,i)));
            col2.Properties.VariableNames{i} = parameters{i};
        end

        %% Merge tables
        TableAll = join(col1,col2);
        %[a,b] = corr(table2array(TableAll),'rows','pairwise');

        %% Get Parameters
        query2 = ['select UKE_Cohort_Parameters.Epi_Nr, T2DnachUKE, Geschlecht_2_w_1_m, OPBMI, Age, OPNASHScore, HbA1c, ISGlucose, ISTriglyceride, Insulin_pmol_L_525_max, HOMA_IR, Cpeptid_pmol_L, Verfettung_Hepatozyten, OPFibrose ',...
            'from UKE_Cohort_Parameters, UKE_Cohort_Scores, UKE_Cohort_FattyContent ',...
            'where UKE_Cohort_Parameters.Epi_Nr = UKE_Cohort_Scores.Epi_Nr and UKE_Cohort_Parameters.Epi_Nr = UKE_Cohort_FattyContent.Epi_Nr'];
        connection = sqlite(dbfile);
        result = fetch(connection,query2);
        close(connection);
        parameters = {'EpiNr', 'T2DnachUKE', 'Geschlecht_2_w_1_m', 'OPBMI', 'Age', 'OPNASHScore','HbA1c','ISGlucose','ISTriglyceride','Insulin','HOMA_IR', 'Cpeptid', 'Hepatolipid','Fibrose'};

        col2 = array2table(zeros(size(result)));
        for i=1:1
            col2(:,i) = array2table(cell2mat(result(:,i)));
            col2.Properties.VariableNames{i} = parameters{i};
        end

        for i=2:14
            col2(:,i) = array2table(str2double(result(:,i)));
            col2.Properties.VariableNames{i} = parameters{i};
        end

        %% Merge tables
        TableAll = join(col1,col2);
        %[a,b] = corr(table2array(TableAll),'rows','pairwise');

        %% Statistics
         % Linear regression
         % T2D w/o NAS
         mdl = fitlm(TableAll(:,[3 4 5 6 2]));
         table_pvalues(k,21) = mdl.Coefficients(2,4);
         table_pvalues(k,22) = mdl.Coefficients(2,1);    
         % BMI 23/24
         mdl = fitlm(TableAll(:,[2 8 4 6 5]));
         table_pvalues(k,23) = mdl.Coefficients(2,4);
         table_pvalues(k,24) = mdl.Coefficients(2,1);     
         % NAS 19/20
         mdl = fitlm(TableAll(:,[7 4 6 5 2]));
         table_pvalues(k,19) = mdl.Coefficients(2,4);
         table_pvalues(k,20) = mdl.Coefficients(2,1);
         % Hepatolipid 27/28
         mdl = fitlm(TableAll(:,[14 4 6 5 2]));
         table_pvalues(k,27) = mdl.Coefficients(2,4);
         table_pvalues(k,28) = mdl.Coefficients(2,1);
         % Glc 15/16
         mdl = fitlm(TableAll(:,[2 4 6 5 9]));
         table_pvalues(k,29) = mdl.Coefficients(2,4);
         table_pvalues(k,30) = mdl.Coefficients(2,1);
         % HbA1c 17/18
         mdl = fitlm(TableAll(:,[2 4 6 5 8]));
         table_pvalues(k,31) = mdl.Coefficients(2,4);
         table_pvalues(k,32) = mdl.Coefficients(2,1);
        % Triglycerides 19/20
         mdl = fitlm(TableAll(:,[2 4 6 5 10]));
         table_pvalues(k,33) = mdl.Coefficients(2,4);
         table_pvalues(k,34) = mdl.Coefficients(2,1);
         % Insulin 21/22
         mdl = fitlm(TableAll(:,[2 4 6 5 11]));
         table_pvalues(k,35) = mdl.Coefficients(2,4);
         table_pvalues(k,36) = mdl.Coefficients(2,1);
         % Cpep 23/24
         mdl = fitlm(TableAll(:,[2 4 6 5 13]));
         table_pvalues(k,37) = mdl.Coefficients(2,4);
         table_pvalues(k,38) = mdl.Coefficients(2,1);
         % HOMA-IR 25/26
         mdl = fitlm(TableAll(:,[12 4 6 5 2]));
         table_pvalues(k,39) = mdl.Coefficients(2,4);
         table_pvalues(k,40) = mdl.Coefficients(2,1);
        
         % Range of Log2 values
         table_pvalues(k,26) = array2table(mean(table2array(TableAll(:,2))));

         %% HUMAN: Is it expressed? 
         connection = sqlite('miRNA.db');
         query3 = ['select hsaObeseLiver_T2D_miRNA4.T2DExpressed ',...
             ' from hsaObeseLiver_T2D_miRNA4, hsa_miRNA4 ',...
             ' where hsaObeseLiver_T2D_miRNA4.ID = hsa_miRNA4.ID ',...
             ' and hsa_miRNA4.TranscriptID like ', term];
         result = fetch(connection,query3);
         close(connection);
         if(result{1}=='F')
             table_pvalues(k,25) = array2table(0);
         else
             table_pvalues(k,25) = array2table(1);
         end
    else
        table_pvalues(k,19:40) = array2table(NaN);
    end
end

table_pvalues = sortrows(table_pvalues,'TranscriptID','ascend');
clear B dev stats mdl Xtable Ytable k i col1 col2 query query2 query3 term term2 termMIR result connection parameters s 
clear Chow CR data Ex H HC HFD i k mdl VSG_d35 VSG_d9 Sham_d35 Sham_d9 term3 dbfile miRNA_IDs TableAll
% Continue with next step
% heatmap