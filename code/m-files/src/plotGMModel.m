function plotGMModel(params, max_muKs, max_rhos, show_labels)
    maxN=length(params.Vs)/2;

    if nargin<4
       show_labels=0; 
    end

    
    

    color_TC=[0     0.56863     0.56863];
    color_WT=[0.83137     0.36863           0];
    
    
    muKs_TC=params.Vs(1:maxN)/max_muKs;
    rhos_TC=params.cs(1:maxN)/max_rhos;
    
    muKs_WT=params.Vs(maxN+1:end)/max_muKs;
    rhos_WT=params.cs(maxN+1:end)/max_rhos;
    
    plot(mean(muKs_TC), mean(rhos_TC), 'o','MarkerFaceColor',color_TC,'Color',color_TC); hold on;
    %ellipse(std(muKs_TC),std(rhos_TC),0,mean(muKs_TC), mean(rhos_TC),color_TC);
    
    plot(mean(muKs_WT), mean(rhos_WT), 'o','MarkerFaceColor',color_WT,'Color',color_WT);
    
    if show_labels==1
        text(mean(muKs_TC), mean(rhos_TC), 'TC','FontSize',20,'HorizontalAlignment','Center','VerticalAlignment','Top','Color',color_TC);
        text(mean(muKs_WT), mean(rhos_WT), 'WT','FontSize',20,'HorizontalAlignment','Center','VerticalAlignment','Bottom','Color',color_WT);
        
    nticks=2;
    else
        nticks=5;
    end
    
    %ellipse(std(muKs_WT),std(rhos_WT),0,mean(muKs_WT), mean(rhos_WT),color_WT);
    
    %xlim([max_muKs/3 max_muKs]);
    %ylim([max_rhos/3 max_rhos]);    
    set(gca,'FontSize',20);


    X_TC = [muKs_TC' rhos_TC'];
    GMModel_TC = fitgmdist(X_TC,1);
    nfitMu_TC=GMModel_TC.mu;
    nfitSigma_TC=GMModel_TC.Sigma;
    gm_TC = gmdistribution(nfitMu_TC,nfitSigma_TC);

    thick = 0.8;   %alpha value
    ch_TC = fcontour(@(x,y)reshape(pdf(gm_TC,[x(:),y(:)]),size(x)),[0 1],'LineWidth',2,'LineColor',color_TC); hold on;
    ch_TC.LevelList=[0 1];
    patches = findobj(ch_TC,'-property','AlphaData');
    for ph = patches
        set(ph,'AlphaData', thick * get(ph,'AlphaData'));
    end
    
    
    
    X_WT = [muKs_WT' rhos_WT'];
    GMModel_WT = fitgmdist(X_WT,1);
    nfitMu_WT=GMModel_WT.mu;
    nfitSigma_WT=GMModel_WT.Sigma;
    gm_WT = gmdistribution(nfitMu_WT,nfitSigma_WT);
    
    ch_WT = fcontour(@(x,y)reshape(pdf(gm_WT,[x(:),y(:)]),size(x)),[0 1],'LineWidth',2,'LineColor',color_WT); hold on;
    ch_WT.LevelList=[0 1];
    patches = findobj(ch_WT,'-property','AlphaData');
    for ph = patches
        set(ph,'AlphaData', thick * get(ph,'AlphaData'));
    end
    


    xticks(0:1/nticks:1);
    xticklabels((0:max_muKs/nticks:max_muKs)/1e-10);

    yticks(0:1/nticks:1);
    yticklabels((0:max_rhos/nticks:max_rhos)/1e8);


    set(gca,'FontSize',20);

    %legend(p, leg,'Location','NorthEastOutside');
    %legend boxoff
    
    %legend([p_WT_K, p_TC_K,p_WT_E, p_TC_E],{'Klebsiella','Klebsiella (plasmid-free)', 'Escherichia', 'Escherichia (plasmid-free)'},'FontSize',16, 'Location','SouthWest');
    %legend boxoff
    xlabel('\mu/K','FontSize',24);
    ylabel('\rho','FontSize',24);