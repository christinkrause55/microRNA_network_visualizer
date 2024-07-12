function [B, dev, stats] = logisticRegressionFromTable(Xtable, Ytable)

    Xreg = table2array(Xtable);
    Yreg = categorical(table2array(Ytable));
    
    [B,dev,stats] = mnrfit(Xreg,Yreg,'interactions','on');
    
    if (size(Xreg,2) == 2)

        t1 = table2array(Ytable) == 1;
        t0 = table2array(Ytable) == 0;
        
        figure
        plot(Xreg(t1,1),Xreg(t1,2),'r*');
        hold on
        plot(Xreg(t0,1),Xreg(t0,2),'bo');
    end
    
end

