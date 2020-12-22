


Mstart=1;
Mend=10;


figure(); 
clf('reset'); 
set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',20); 


M_invasion=zeros(length(Ms),length(Cs));
M_pf=zeros(length(Ms),length(Cs));

minColor=0.;
maxColor=1;


color_TC=[0     0.56863     0.56863];
color_WT=[0.83137     0.36863           0];
cmap_freq=makeColorMap(color_WT,[1 1 1],color_TC, 201);
colormap(cmap_freq);

for sni=Mstart:length(Ms)

    %percent_unstable=zeros(1,length(Ms));
    %mean_unstable=zeros(1,length(Ms));
    
    for ci=1:length(Cs)
        %percent_unstable(s)=sum(all_t_end(s,:, c)==-1)/numExperiments;
        %mean_unstable(s)=mean(all_t_end(all_t_end(s,:, c)>-1));
        
        invaded=0;
        pf=0;
        for ni=1:numExperiments
           %success_invasnion 
           if all_t_end(sni, ci, ni)==-1
              invaded=invaded+1; 
           end
           pf=pf+all_pf(sni, ci, ni)/numExperiments;
        end
        invasion_rate=invaded/numExperiments;
        %icolor=floor(invasion_rate*(numColors-1))+1;
        %this_color=cmap(icolor,:);
        
        M_invasion(sni, ci+1)=invasion_rate;
        M_pf(sni, ci+1)=pf;
        
        if pf<minColor
            this_color=[1 1 1];
        else
        
            pfColor=(pf-minColor)/(maxColor-minColor);
            if pfColor<0
                pfColor=0;
            elseif pfColor>1
                pfColor=1; 
            end

            icolor=floor(pfColor*(numColors-1))+1;
            this_color=cmap_freq(icolor,:);
        end
           
        padd=.0;
        rectangle('Position',[sni-.5+padd,ci-.5+padd, 1-padd, 1-padd],'FaceColor',this_color,'Curvature',[0 0]); hold on;
        
    end

    % PLOT STABILITY
    %plot(1:length(Ms), (1-percent_unstable),'o:'); hold all;
   
end
%legend(Cs);
axis tight
    xticks(Mstart:1:Mend);
    xticklabels(Mstart:1:Mend);
    xlim([Mstart-.5, Mend+.5]);

yticks(1:2:length(Cs));
str_yticks=['0',strcat(strcat('10^{',strsplit(num2str(log10(Cs(Cs>0))),' ')),'}')];
yticklabels(str_yticks);

colormap(cmap_freq);
cbh=colorbar();
yt=linspace(0, 1, 3);
yt_lbl=num2cell(linspace(minColor, maxColor, 3)*100);

if maxColor<1
    yt_lbl{end}=['>',num2str(yt_lbl{end})]; 
end

if minColor>0 
    yt_lbl{1}=['<',num2str(yt_lbl{end})];
    
end

set(cbh,'YTick',yt);
set(cbh,'YTickLabel',yt_lbl);
ylabel(cbh,'Plasmid-bearing cells (%)','FontSize',24);
%xlim([0.4 length(Ms)+.6]);
%ylim([0.4 length(Cs)+.6]);
%ylim([0, 1]);
xlabel('Number of strains in community','FontSize',24);
ylabel('Conjugation rate','FontSize',24);
set(gca,'FontSize',20)
