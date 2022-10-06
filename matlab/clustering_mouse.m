%% This skript needs TraitTable2 and 3 generated for Clustering
% Generate data. 
[coeff,score,latent,tsquared,explained,mu] = pca(table2array(traitTable(:,1:end-1)));

figure(1)
for k =1:9
    % kMeans
    idxSpectral = kmeans(table2array(traitTable(:,1:end-1)),k);
    group = zeros(54,k+9);
    for i=1:k
        group(:,i) = idxSpectral == i;
    end

    colors = ['r'; 'b'; 'g'; 'm'; 'k'; 'y'; 'c'; 'r'; 'b'; 'g'];

    %Each Group
    group(:,end-8) = strcmp(string(table2array(traitTable(:,end))),"chow");
    group(:,end-7) = strcmp(string(table2array(traitTable(:,end))),"HFD");
    group(:,end-6) = strcmp(string(table2array(traitTable(:,end))),"H>C");
    group(:,end-5) = strcmp(string(table2array(traitTable(:,end))),"CR");
    group(:,end-4) = strcmp(string(table2array(traitTable(:,end))),"Ex4");
    group(:,end-3) = strcmp(string(table2array(traitTable(:,end))),"Sham d9");
    group(:,end-2) = strcmp(string(table2array(traitTable(:,end))),"VSG d9");
    group(:,end-1) = strcmp(string(table2array(traitTable(:,end))),"Sham d35");
    group(:,end-0) = strcmp(string(table2array(traitTable(:,end))),"VSG d35");

    subplot(3,3,k)
    for i=1:k
        t = group(:,i) == 1;
        scatter(score(t,1),score(t,2),append(colors(i),'o'))
        chow = group(:,i) == 1 & group(:,end-8) == 1;
        hfd = group(:,i) == 1 & group(:,end-7) == 1;
        hc = group(:,i) == 1 & group(:,end-6) == 1;
        cr = group(:,i) == 1 & group(:,end-5) == 1;
        ex = group(:,i) == 1 & group(:,end-4) == 1;
        sd9 = group(:,i) == 1 & group(:,end-3) == 1;
        vd9 = group(:,i) == 1 & group(:,end-2) == 1;
        sd35 = group(:,i) == 1 & group(:,end-1) == 1;
        vd35 = group(:,i) == 1 & group(:,end-0) == 1;

        hold on
        scatter(score(chow,1),score(chow,2),'filled',colors(i))
    end
    title("kmeans clustering with k = " + k)
    xlabel("PC1 (" + string(explained(1)) + " %)")
    ylabel("PC2 (" + string(explained(2)) + " %)")
    hold off
end


%% Clear workspace
clear col2 colors connection dbfile EpiNr i k latent mu nd orderND orderT2D 
clear parameters query query2 result s score FCup selectFCupdw 
clear t t2 t2d take term term2 tmp toi values_twist traitTable2 traitTable3 selectFCup explained
clear test  mean_values_clustering tsquared coeff idxSpectral group

% drawGraph