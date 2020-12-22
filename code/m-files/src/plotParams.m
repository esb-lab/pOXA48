function plotParams(params, max_muKs, max_rhos, nfig)

    %if nargin<4
    %   nfig=2; 
    %end
    %figure(nfig); %clf('reset');
    %set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white');

    nparams_muKs=params.Vs./max_muKs;
    nparams_rhos=params.cs./max_rhos;

    p=[];
    leg={};
    for istrain=1:params.numStrains

%        if ~strcmp(params.strains{istrain},'X')

            if isfield(params,'colors')
                this_color=params.colors(istrain,:);
            else
                this_color=[rand rand rand];%[.7 .7 .7];
            end
            
            if isfield(params,'species')
                %Klebsiella
                if strcmp(params.species(istrain),'K')
                    str_marker='d';
                else %E coli
                    str_marker='o';
                end
            else
                str_marker='o';
            end
            plot(nparams_muKs(istrain+params.numStrains),nparams_rhos(istrain+params.numStrains), str_marker,'LineWidth',1,'MarkerEdgeColor',this_color,'MarkerSize',10); hold on;
            p(istrain)=plot(nparams_muKs(istrain), nparams_rhos(istrain), str_marker,'MarkerFaceColor',this_color,'LineWidth',1,'MarkerEdgeColor',this_color,'MarkerSize',10); hold on;
            
            if isfield(params,'strains')
                leg{istrain}=extractBefore(params.strains{istrain},"_");
            end
 %       end
    end
    
    
    nticks=5;
    axis([0.2 1 0.4 1]);

    xticks(0:1/nticks:1);
    xticklabels(0:max_muKs/nticks:max_muKs);

    yticks(0:1/nticks:1);
    yticklabels(0:max_rhos/nticks:max_rhos);


    set(gca,'FontSize',20);

    %legend(p, leg,'Location','NorthEastOutside');
    %legend boxoff
    
    %legend([p_WT_K, p_TC_K,p_WT_E, p_TC_E],{'Klebsiella','Klebsiella (plasmid-free)', 'Escherichia', 'Escherichia (plasmid-free)'},'FontSize',20, 'Location','SouthWest');
    %legend boxoff
    xlabel('Specific affinity (V_{max}/K_m)','FontSize',24);
    ylabel('Cell efficiency (\rho)','FontSize',24);