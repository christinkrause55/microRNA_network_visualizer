%% The following lines need a p-Value table generated by callDatabase.m
%% Setup
% this database is availabe upon request
dbfile = 'GeneExpressions.db';
s='''';
%% Select miRNA with association of interest
%

select = table2array(table_pvalues(:,5)) == 1 & table2array(table_pvalues(:,6)) > 2.3 & ...
        (table2array(table_pvalues(:,27)) < 0.05 | table2array(table_pvalues(:,27)) < 0.05);    
        %table2array(table_pvalues(:,15)) < 0.05;   %Glucose
        %table2array(table_pvalues(:,17)) < 0.05;   %HbA1c
        %table2array(table_pvalues(:,19)) < 0.05;   %Triglycerides
        %table2array(table_pvalues(:,21)) < 0.05;   %Insulin
        %table2array(table_pvalues(:,23)) < 0.05;   %Cpeptide
        %table2array(table_pvalues(:,25)) < 0.05;   %HOMA-IR
        %table2array(table_pvalues(:,27)) < 0.05;   %NAS
        %table2array(table_pvalues(:,7)) < 0.05;    %T2D

disp(sum(select) + " miRNAs to process")

toi = table_pvalues(select,:);

% Exclude precursor mir
take = ~contains(table2array(toi(:,end)),'r');
tmp = toi(take,:);
toi = tmp;

% Setup tables
values = zeros(size(table2array(toi(:,4)),1),40);
IDs = cell(size(table2array(toi(:,4)),1),1);
EpiNr = zeros(1,40);

% Fill tables by miRNA transcript ID
for k=1:size(table2array(toi(:,4)),1)
    connection = sqlite(dbfile);
    tmp = table2array(toi(k,end));
    term = strcat(s,tmp{1},s);
    term2 = strcat(s,'sample',s);
    
    query = ['select * from values_miR4_col where values_miR4_col.TranscriptID like ', ...
        term, ' or values_miR4_col.TranscriptID = ',...
        term2];
    result = fetch(connection,query);
    close(connection);
    values(k,:) = cell2mat(result(2,3:end));
    IDs(k,1) = result(2,2);
    EpiNr = cell2mat(result(1,3:end));
end

%% Merge miRNA values with phenotypes of interest
tmp = [EpiNr; values];
values_twist = array2table(tmp.');
values_twist.Properties.VariableNames{1} = 'EpiNr';

for i=2: size(table2array(values_twist(2,:)),2)
    values_twist.Properties.VariableNames{i} = IDs{i-1};
end

query2 = ['select ID, T2D from Cohort_Parameters'];
    connection = sqlite(dbfile);
    result = fetch(connection,query2);
    close(connection);
    parameters = {'ID', 'T2D'};

    col2 = array2table(zeros(size(result)));
    for i=1
        col2(:,i) = array2table(cell2mat(result(:,i)));
        col2.Properties.VariableNames{i} = parameters{i};
    end
    for i=2:size(result,2)
        col2(:,i) = array2table(str2double(result(:,i)));
        col2.Properties.VariableNames{i} = parameters{i};
    end

%% Merge tables and draw Heatmap
traitTable = join(values_twist,col2);
traitTable = sortrows(traitTable,'T2D','ascend');

% Order for ND
clustergND = clustergram(table2array(traitTable(1:20,2:(end-1))),'Colormap',redbluecmap,...
'Standardize','Column', 'RowLabels',table2array(traitTable(1:20,1)),'ColumnLabels',IDs, ...
'Linkage','complete','Dendrogram',3);
orderND = array2table(str2num(cell2mat(clustergND.RowLabels)));
orderND.Properties.VariableNames{1} = 'ID';
tempTable1 = join(orderND,traitTable(1:20,:));

% Order for T2D
clustergT2D = clustergram(table2array(traitTable(21:40,2:(end-1))),'Colormap',redbluecmap,...
'Standardize','Column', 'RowLabels',table2array(traitTable(21:40,1)),'ColumnLabels',IDs, ...
'Linkage','complete','Dendrogram',3);
orderT2D = array2table(str2num(cell2mat(clustergT2D.RowLabels)));
orderT2D.Properties.VariableNames{1} = 'ID';
tempTable2 = join(orderT2D,traitTable(21:40,:));
% Merge Tables
traitTable = [tempTable1; tempTable2];

clusterg = clustergram(table2array(traitTable(:,2:(end-1))),'Colormap',redbluecmap,...
    'Standardize','Column', 'RowLabels',table2array(traitTable(:,end)),'ColumnLabels',IDs, ...
    'Linkage','complete','Dendrogram',4,'Cluster','row');
set(0,'ShowHiddenHandles','on')
allhnds = get(0,'Children');
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 16)
cgfig = findall(0,'type','figure', 'tag', 'Clustergram'); % Figure handle
cgax = findobj(cgfig, 'type','axes','tag','HeatMapAxes'); % main axis handle
%cgax.Fontsize = 20;

%% Sort by FoldChange
% shorten Table to exclude ID and group
traitTable2 = traitTable(:,2:end);
t2 = table2array(traitTable2(:,end)) == 1;
nd = table2array(traitTable2(:,end)) == 0;
traitTable3 = traitTable(:,2:end-1);

FC = zeros(size(traitTable3(1,:),2),4);
for i=1:size(FC,1)
    FC(i,1) = tukeyBiweight2(table2array(traitTable2(t2,i)),5,0.0001)/tukeyBiweight2(table2array(traitTable2(nd,i)),5,0.0001);
    %FC(i,1) = mean(table2array(traitTable3(t2,i)))/mean(table2array(traitTable3(nd,i)));
    FC(i,2) = (std(table2array(traitTable3(t2,i)))/sqrt(20))/tukeyBiweight2(table2array(traitTable2(nd,i)),5,0.0001);
    FC(i,3) = 1;
    FC(i,4) = (std(table2array(traitTable3(nd,i)))/sqrt(20))/tukeyBiweight2(table2array(traitTable2(nd,i)),5,0.0001);
end

clear t2 nd orderT2D dbfile s specialTargets toi tmp take values EpiNr cgfig cgax ans coeff
clear clusterg clusterND clusterT2D col2 i k orderND parameters query query2 result
clear selectFCup selectFCupdw term term2 tempTable1 tempTable2 valuesTwist
clear clustergND clustergT2D 
clustering
 

    