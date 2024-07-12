%% This skript needs TraitTable2 and 3 generated for Clustering
% Figures from this script were not used in publication
% Generate data. 
[coeff,score,latent,tsquared,explained,mu] = pca(table2array(traitTable3));

figure(1)
title("Spectral clustering")
for k =1:6
    % Spectral
    idxSpectral = spectralcluster(table2array(traitTable3),k);
    group = zeros(40,k+2);
    for i=1:k
        group(:,i) = idxSpectral == i;
    end

    colors = ['r'; 'b'; 'g'; 'm'; 'k'; 'y'];

    %T2D individuals
    group(:,end-1) = table2array(traitTable2(:,end)) == 1;
    %ND individuals
    group(:,end) = table2array(traitTable2(:,end)) == 0;

    subplot(3,2,k)
    for i=1:k
    t = group(:,i) == 1;
    scatter(score(t,1),score(t,2),append(colors(i),'o'))
    t2d = group(:,i) == 1 & group(:,end-1) == 1;
    nd = group(:,i) == 1 & group(:,end) == 1;
    hold on
    scatter(score(nd,1),score(nd,2),'filled',colors(i))
    end
    title("Spectral clustering with k = " + k)
    xlabel("PC1")
    ylabel("PC2")
    hold off
end

figure(2)
for k =1:6
    % kMeans
    idxSpectral = kmeans(table2array(traitTable3),k);
    group = zeros(40,k+2);
    for i=1:k
        group(:,i) = idxSpectral == i;
    end

    colors = ['r'; 'b'; 'g'; 'y'; 'm'; 'k'];

    %T2D individuals
    group(:,end-1) = table2array(traitTable2(:,end)) == 1;
    %ND individuals
    group(:,end) = table2array(traitTable2(:,end)) == 0;

    subplot(3,2,k)
    for i=1:k
    t = group(:,i) == 1;
    scatter(score(t,1),score(t,2),append(colors(i),'o'))
    t2d = group(:,i) == 1 & group(:,end-1) == 1;
    nd = group(:,i) == 1 & group(:,end) == 1;
    hold on
    scatter(score(nd,1),score(nd,2),'filled',colors(i))
    end
    title("kmeans clustering with k = " + k)
    xlabel("PC1")
    ylabel("PC2")
    hold off
end

% kMeans sep
k=3;
    idxSpectral = kmeans(table2array(traitTable3),k);
    group = zeros(40,k+2);
    for i=1:k
        group(:,i) = idxSpectral == i;
    end

    colors = ['r'; 'b'; 'g'; 'y'; 'm'; 'k'];

    %T2D individuals
    group(:,end-1) = table2array(traitTable2(:,end)) == 1;
    %ND individuals
    group(:,end) = table2array(traitTable2(:,end)) == 0;

figure
mean_values_clustering = zeros(2,k);
    for i=1:k
        t = group(:,i) == 1;
        scatter(score(t,1),score(t,2),append(colors(i),'o'))
        t2d = group(:,i) == 1 & group(:,end-1) == 1;
        nd = group(:,i) == 1 & group(:,end) == 1;
        mean_values_clustering(1,i) = sum(t2d);
        mean_values_clustering(2,i) = sum(nd);
        hold on
        scatter(score(nd,1),score(nd,2),'filled',colors(i))
    end
    title("kmeans clustering with k = " + k)
    xlabel("PC1 (" + string(explained(1)) + " %)")
    ylabel("PC2 (" + string(explained(2)) + " %)")
    hold off

    mean_values_clustering

%% Clear workspace
clear col2 colors connection dbfile EpiNr i k latent mu nd orderND orderT2D 
clear parameters query query2 result s score FCup selectFCupdw 
clear t t2 t2d take term term2 tmp toi values_twist traitTable2 traitTable3 selectFCup explained
clear test  mean_values_clustering tsquared coeff idxSpectral group

drawGraph
