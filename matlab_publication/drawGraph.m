%% Set miRNA of interest
values_T2D = values_NASH;


mRNAs = {'SCD','FASN','ELOVL6','ACACA','PNPLA2','LDLR'...
    'PPARA','LRP6',...
    'PDK1','INSR','IRS1','IRS2','FOXO1','PCK1','G6PC','SLC2A2'};

% Take only Trait-Associated miRNAs generated from heatmap and which
% where clustered by clustering
miRNA_List = IDs;
miRNAs = [];
for i=1:size(miRNA_List,1)
    indx = find(strcmp(table2array(values_T2D(:,3)),miRNA_List{i}));
    if(size(indx,1) ~= 0 )
        % Get positions of miRNAs in adjacency/miRNA List and save them
        % in miRNAs
        miRNAs = [miRNAs indx(1,1)];
    else
        disp("not found " + miRNA_List{i})
        miRNAs = [miRNAs NaN];
    end
end

% Needs miRNAs, target genes, miRNA sequences, mRNA sequences
matches = callSeedPrediction(table2array(values_T2D(miRNAs.',3)),mRNAs,human_miR, human_mRNA);
GOI = [];
  % get Information about possible selected target genes
  for i=1:size(mRNAs,2)
    indx = find(strcmp(table2array(values_Gene(:,3)),mRNAs{i}));
    if(size(indx,1) ~= 0 )
    GOI = [GOI indx(1,1)];
  end
end
adjacency_select = adjacency_matrix(miRNAs,GOI);


%% Add prediction values to adjacency_select
pred = matches ~= 0;
adjacency_select(pred) = 4;

% Delete those miRNAs without any association
miRNAs_t = [];
for i=1:size(miRNAs,2)
    disp(i)
    if (~(sum(adjacency_select(i,:) ~= 0) == 0))
        miRNAs_t = [miRNAs_t i];
    end
    
end
adjacency_select = adjacency_select(miRNAs_t,:);

miRNA_Expression = traitTable(:,[1 miRNAs_t+1]);
mRNA_Expression = getGraphWeights(traitTable,mRNAs);

adjacency_pValues = zeros(size(adjacency_select));
adjacency_rValues = zeros(size(adjacency_select));

% Correlate miRNA expression with gene expression
for row=1:size(adjacency_select,1)
  for col=1:size(adjacency_select,2)
        [adjacency_rValues(row,col),adjacency_pValues(row,col)] = corr(table2array(miRNA_Expression(:,row+1)),table2array(mRNA_Expression(:,col+1)),'rows','pairwise');
        adjacency_rValues(row,col) = adjacency_rValues(row,col)*-1;
  end
end

% Cleanup adjacency_select for non used genes
 tmp = [];
 tmp2 = [];
 tmp3 = [];
 GOI_tmp = [];
 for i=1:size(adjacency_select,2)
     if(sum(adjacency_select(:,i)) ~= 0)
         tmp = [tmp adjacency_select(:,i)];
         tmp2 = [tmp2 adjacency_pValues(:,i)];
         tmp3 = [tmp3 adjacency_rValues(:,i)];
         GOI_tmp = [GOI_tmp GOI(i)];
     end
 end
 
adjacency_select = tmp;
adjacency_matrix_p = tmp2;
adjacency_rValues = tmp3;
GOI = GOI_tmp;
 
nLabels = [string(IDs(miRNAs_t,1)); string(table2cell(values_Gene(GOI,3)))].';

%% plot
negSource = [];
negTarget = [];
nSource = [];
nTarget = [];
weight = [];
eLabels = [];
nLabels = [string(IDs(miRNAs_t,1)); string(table2cell(values_Gene(GOI,3)))].';

for i=1:size(adjacency_select,1)
    if(sum(adjacency_select(i,:)) == 0)
            continue
    end
    for j=1:size(adjacency_select,2)
        if(adjacency_select(i,j) ~= 0)
            nSource = [nSource i];
            nTarget = [nTarget j+size(adjacency_select,1)];
            eLabels = [eLabels string(adjacency_rValues(i,j))];  
            if(adjacency_rValues(i,j) < 0)
                negSource = [negSource i];
                negTarget = [negTarget j+size(adjacency_select,1)];
                weight = [weight abs(adjacency_rValues(i,j))*10];
            else
                weight = [weight 0.1];
            end
        end
    end
end

G2 = graph(nSource, nTarget);
p = plot(G2,'Layout','layered','NodeFontSize',12,'LineWidth',weight,'NodeLabel',nLabels,'EdgeColor',[0.1 0.1 0.1])
p.Marker = 'o';
p.NodeColor = 'r';
p.MarkerSize = 8;
highlight(p,[size(adjacency_select,1)+1:1:size(adjacency_select,1)+size(adjacency_select,2)],'NodeColor','k','Marker', 's','MarkerSize',10)
highlight(p,negSource,negTarget,'EdgeColor',[0.25 0.5 1])

clear indx isNegative indxt tmpWidth i tmp row col GOI miRNA_t w miRNA_List
clear RNAseq_mode eWidth eLabels k last RNAseq_targeted
clear nSource nTarget nLabels eWidth indxt eLabels firstFlag miRNAs_t G p
effectSizes
