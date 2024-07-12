% Gene expression within this database is available upon request
dbfile = 'GeneExpressions.db';

%% Get miRNA Names
s='''';
termMIR = strcat(s,'hsa-%',s);
connection = sqlite(dbfile);
query = ['select ID,TranscriptID from values_miR4_col where values_miR4_col.TranscriptID like ', ...
    termMIR];
miRNA_IDs = fetch(connection,query);
close(connection);

table_pvalues = array2table(zeros(size(miRNA_IDs,1),34));
table_pvalues.Properties.VariableNames{1} = 'pT2D';
table_pvalues.Properties.VariableNames{2} = 'eT2D';
table_pvalues.Properties.VariableNames{3} = 'pTSD_NAS';
table_pvalues.Properties.VariableNames{4} = 'eT2D_NAS';
table_pvalues.Properties.VariableNames{5} = 'isexpressed';
table_pvalues.Properties.VariableNames{6} = 'meanAll';
table_pvalues.Properties.VariableNames{7} = 'p_logistic';
table_pvalues.Properties.VariableNames{8} = 'e_logistic';
table_pvalues.Properties.VariableNames{9} = 'p_logistic_NAS';
table_pvalues.Properties.VariableNames{10} = 'e_logistic_NAS';
table_pvalues.Properties.VariableNames{11} = 'pAge';
table_pvalues.Properties.VariableNames{12} = 'eAge';
table_pvalues.Properties.VariableNames{13} = 'pBMI';
table_pvalues.Properties.VariableNames{14} = 'eBMI';
table_pvalues.Properties.VariableNames{15} = 'pGlc';
table_pvalues.Properties.VariableNames{16} = 'eGlc';
table_pvalues.Properties.VariableNames{17} = 'pHbA1c';
table_pvalues.Properties.VariableNames{18} = 'eHbA1c';
table_pvalues.Properties.VariableNames{19} = 'pTriglycerides';
table_pvalues.Properties.VariableNames{20} = 'eTriglycerides';
table_pvalues.Properties.VariableNames{21} = 'pInsulin';
table_pvalues.Properties.VariableNames{22} = 'eInsulin';
table_pvalues.Properties.VariableNames{23} = 'pCpep';
table_pvalues.Properties.VariableNames{24} = 'eCpep';
table_pvalues.Properties.VariableNames{25} = 'pHOMAIR';
table_pvalues.Properties.VariableNames{26} = 'eHOMAIR';
table_pvalues.Properties.VariableNames{27} = 'pNAS';
table_pvalues.Properties.VariableNames{28} = 'eNAS';
table_pvalues.Properties.VariableNames{29} = 'pHepatolipid';
table_pvalues.Properties.VariableNames{30} = 'eHepatolipid';
table_pvalues.Properties.VariableNames{31} = 'pFibrose';
table_pvalues.Properties.VariableNames{32} = 'eFibrose';
table_pvalues.Properties.VariableNames{33} = 'pSex';
table_pvalues.Properties.VariableNames{34} = 'eSex';

table_pvalues(:,35) = array2table(miRNA_IDs(:,2));
table_pvalues.Properties.VariableNames{35} = 'TranscriptID';
table_pvalues(:,36) = array2table(miRNA_IDs(:,1));
table_pvalues.Properties.VariableNames{36} = 'ArrayID';

for k=1:size(miRNA_IDs,1)
    disp('moving to position ' + string(k))
    
    %% Get miRNA Information
    connection = sqlite(dbfile);
    term = strcat(s,miRNA_IDs{k},s);
    term2 = strcat(s,'sample',s);

    query = ['select * from values_miR4_col where values_miR4_col.ID like ', ...
        term, ' or values_miR4_col.TranscriptID = ',...
        term2];

    result = fetch(connection,query);
    close(connection);

    col1 = array2table(cell2mat(result(:,3:end)).');
    col1.Properties.VariableNames{1} = 'EpiNr';
    col1.Properties.VariableNames{2} = miRNA_IDs{k,2};

    %% Get Parameters
    query2 = ['select Cohort_Parameters.ID, T2D, sex, BMI, age, NASscore, HbA1c, glucose, triglyceride, insulin, HOMA_IR, Cpeptid, fat, fibrosis ',...
        'from Cohort_Parameters];
    connection = sqlite(dbfile);
    result = fetch(connection,query2);
    close(connection);
    parameters = {'ID', 'T2D', 'sex', 'BMI', 'age', 'NASscore','HbA1c','glucose','triglyceride','insulin','HOMA_IR', 'Cpeptid', 'fat','fibrosis'};

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
     table_pvalues(k,1) = mdl.Coefficients(2,4);
     table_pvalues(k,2) = mdl.Coefficients(2,1);
     
     % T2D w/ NAS
     mdl = fitlm(TableAll(:,[3 4 5 6 7 2]));
     table_pvalues(k,3) = mdl.Coefficients(2,4);
     table_pvalues(k,4) = mdl.Coefficients(2,1);
     
     % Age 11/12
     mdl = fitlm(TableAll(:,[2 8 4 5 7 6]));
     table_pvalues(k,11) = mdl.Coefficients(2,4);
     table_pvalues(k,12) = mdl.Coefficients(2,1);
     
     % BMI 13/14
     mdl = fitlm(TableAll(:,[2 8 4 6 7 5]));
     table_pvalues(k,13) = mdl.Coefficients(2,4);
     table_pvalues(k,14) = mdl.Coefficients(2,1);
     
     % Glc 15/16
     mdl = fitlm(TableAll(:,[2 4 6 5 9]));
     table_pvalues(k,15) = mdl.Coefficients(2,4);
     table_pvalues(k,16) = mdl.Coefficients(2,1);
     
     % HbA1c 17/18
     mdl = fitlm(TableAll(:,[2 4 6 5 7 8]));
     table_pvalues(k,17) = mdl.Coefficients(2,4);
     table_pvalues(k,18) = mdl.Coefficients(2,1);
     
     % Triglycerides 19/20
     mdl = fitlm(TableAll(:,[2 4 6 5 10]));
     table_pvalues(k,19) = mdl.Coefficients(2,4);
     table_pvalues(k,20) = mdl.Coefficients(2,1);
     
     % Insulin 21/22
     mdl = fitlm(TableAll(:,[2 4 6 5 11]));
     table_pvalues(k,21) = mdl.Coefficients(2,4);
     table_pvalues(k,22) = mdl.Coefficients(2,1);
     
     % Cpep 23/24
     mdl = fitlm(TableAll(:,[2 4 6 5 13]));
     table_pvalues(k,23) = mdl.Coefficients(2,4);
     table_pvalues(k,24) = mdl.Coefficients(2,1);
     
     % HOMA-IR 25/26
     mdl = fitlm(TableAll(:,[12 4 6 5 2]));
     table_pvalues(k,25) = mdl.Coefficients(2,4);
     table_pvalues(k,26) = mdl.Coefficients(2,1);
     
     % NAS 27/28
     mdl = fitlm(TableAll(:,[7 4 6 5 2]));
     table_pvalues(k,27) = mdl.Coefficients(2,4);
     table_pvalues(k,28) = mdl.Coefficients(2,1);
     
     % Hepatolipid 29/30
     mdl = fitlm(TableAll(:,[14 4 6 5 2]));
     table_pvalues(k,29) = mdl.Coefficients(2,4);
     table_pvalues(k,30) = mdl.Coefficients(2,1);
     
     % Fibrosis 31/32
     mdl = fitlm(TableAll(:,[15 4 6 5 2]));
     table_pvalues(k,31) = mdl.Coefficients(2,4);
     table_pvalues(k,32) = mdl.Coefficients(2,1);
     
     % Sex 33/34
     mdl = fitlm(TableAll(:,[2 8 6 5 7 4]));  
     Xtable = TableAll(:,[2 8 6 5 7]);
     Ytable = TableAll(:,4);
     [B, dev, stats] = logisticRegressionFromTable(Xtable, Ytable);
     
     table_pvalues(k,33) = array2table(stats.p(2,:));
     table_pvalues(k,34) = array2table(stats.beta(2,:));

     % Range of Log2 values
     table_pvalues(k,6) = array2table(mean(table2array(TableAll(:,2))));
     
     % Logistic regression analysis
     Xtable = TableAll(:,[2 4 5 6]);
     Ytable = TableAll(:,3);
     [B, dev, stats] = logisticRegressionFromTable(Xtable, Ytable);
     
     table_pvalues(k,7) = array2table(stats.p(2,:));
     table_pvalues(k,8) = array2table(stats.beta(2,:));
     
     Xtable = TableAll(:,[2 4 5 6 7]);
     Ytable = TableAll(:,3);
     [B, dev, stats] = logisticRegressionFromTable(Xtable, Ytable);
     
     table_pvalues(k,9) = array2table(stats.p(2,:));
     table_pvalues(k,10) = array2table(stats.beta(2,:));
     
     % Is it expressed?
     connection = sqlite('miRNA.db');
     query3 = ['select hsaObeseLiver_T2D_miRNA4.T2DExpressed ',...
         ' from hsaObeseLiver_T2D_miRNA4, hsa_miRNA4 ',...
         ' where hsaObeseLiver_T2D_miRNA4.ID = hsa_miRNA4.ID ',...
         ' and hsa_miRNA4.ID like ', term];
     result = fetch(connection,query3);
     close(connection);
     if(result{1}=='F')
         table_pvalues(k,5) = array2table(0);
     else
         table_pvalues(k,5) = array2table(1);
     end
end

clear B dev stats mdl Xtable Ytable k i col1 col2 query query2 query3 term term2 termMIR result connection parameters s 
% Continue with next step
% heatmap
