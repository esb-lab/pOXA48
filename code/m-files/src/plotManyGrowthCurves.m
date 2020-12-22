function plotManyGrowthCurves(time, y, params, max_muKs, max_rhos, gm, dyn)


if nargin<3
    strain_labels={};
    strain_colors=jet(params.numStrains);
else
    if isfield(params,'strains')
        strain_labels=params.strains;
    else
        strain_labels={};
    end
    
    if isfield(params,'colors')
        strain_colors=params.colors;
    else
       strain_colors=jet(params.numStrains);
    end
end


if nargin<4
    max_muKs=10e-10;
    max_rhos=12e8;
end


if nargin<6
    gm=[];
end

if nargin<7
   dyn=nan; 
end

%level_set=[1];
level_color=[.8 .8 .8];
rmax=2400;
rmin=24;


figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'fontsize',20);
set(gcf, 'Units','normalized','Position',[0. 0.2 0.6 .3]);

T=time(end);

S=y(:,1);
B=y(:,2:end);


[~, totStrains]=size(B);
numStrains=totStrains/2;

tot=sum(B,2);  %Total density of bacteria

maxY=max([max(sum(B(:,1:numStrains),2)), max(sum(B(:,numStrains+1:end),2)) ]);
labels={};
%p_prop=zeros(1,numStrains);

% Plot contour of metacommunity
subaxis(1,3,1,'SpacingHoriz',0.075,'PaddingBottom',.1);
if ~isempty(gm)
    fcontour(@(x,y)reshape(pdf(gm,[x(:),y(:)]),size(x)),[0 1],'LineWidth',.5,'LineColor',level_color,'MeshDensity',100); hold on;
end

for i=1:numStrains
    
%    if ~strcmp(strain_labels{i},'X')
        
        if tot(end)>0

            this_color=strain_colors(i,:);

            %Plot parameters w/ circle radius proportional to relative density
            subaxis(1,3,1,'SpacingHoriz',0.075,'PaddingBottom',.1);
            
            this_freqTC=B(end,i)./tot(end);
            this_freqWT=B(end,i+params.numStrains)./tot(end);

            %Klebsiella
            %if strcmp(params.species(i),'K')
            %    str_marker='d';
            %else %E coli
                str_marker='o';
            %end

            rmarker_TC=rmin+(rmax-rmin)*this_freqTC;
            rmarker_WT=rmin+(rmax-rmin)*this_freqWT;
            
            
            %Plasmid-bearing
            %scatter(params.Vs(i)/max_muKs, params.cs(i)/max_rhos,rmarker_WT,this_color,'filled'); hold on;
            s_TC=scatter(params.Vs(i)/max_muKs, params.cs(i)/max_rhos,rmarker_TC,this_color); hold on;
            s_TC.LineWidth = 1;
            s_TC.MarkerEdgeColor = this_color;
            s_TC.MarkerFaceColor = this_color;
            s_TC.MarkerFaceAlpha = 0.5;
            
            %Plasmid-free
            s_WT=scatter(params.Vs(i+numStrains)/max_muKs, params.cs(i+numStrains)/max_rhos, rmarker_WT,this_color,'MarkerFaceAlpha',1); hold on;
            s_WT.LineWidth = 2;
            s_WT.MarkerEdgeColor = this_color;
            s_WT.MarkerFaceColor = [1 1 1];
            s_WT.MarkerFaceAlpha = 0.01;
            

            %if B(1,i)>0  %Plasmid-free at start of experiment
            %    plot(params.Vs(i)/max_muKs, params.cs(i)/max_rhos, '.k','MarkerFaceColor','k','MarkerSize',8); hold on;
            %end
                
            %plot(params.Vs(i+params.numStrains)/max_muKs,params.cs(i+params.numStrains)/max_rhos, str_marker,'LineWidth',2,'MarkerEdgeColor',this_color,'MarkerSize',rmarker_WT); hold on;
            %labels{i}=extractBefore(strain_labels{i},'_');

            set(gca,'fontsize',20);
            title(['M=',num2str(params.numStrains)],'FontSize',20);
            xlabel('Specific affinity (\mu/K)','FontSize',24);
            ylabel('Cell efficiency (\rho)','FontSize',24);
            nticks=5;
            axis([0 1 0 1]);
            xticks(0:1/nticks:1);
            xticklabels(0:max_muKs/nticks:max_muKs);
            yticks(0:1/nticks:1);
            yticklabels(0:max_rhos/nticks:max_rhos);

            %Plot density CT vs WT
            subaxis(1,3,2,'SpacingHoriz',0.075);
            %semilogy(time, B(:,i),'-', 'Color', strain_colors(i,:)); hold on;
            %semilogy(time, B(:,i+numStrains),'--', 'Color', strain_colors(i,:)); hold on;
            %yticks(10.^linspace(0,9,10));
            
            plot(time, B(:,i),'-', 'Color', strain_colors(i,:)); hold on;
            plot(time, B(:,i+numStrains),'--', 'Color', strain_colors(i,:)); hold on;
            ylim([0 10^9]);
            
            if ~isnan(dyn)
                 if dyn==-1
                     title(['Plasmid extinction'],'FontSize',20);
                 elseif dyn==0
                     title(['Plasmid co-existence'],'FontSize',20);
                 elseif dyn==1
                     title(['Plasmid fixation'],'FontSize',20);
                 end
            end
            
            
            ylabel('Density','fontsize',24);
            axis([0 time(end) 1 2*maxY]);
            xlabel('Time','fontsize',24);
            if T==24
                xticks(0:T/6:T);
            end
            set(gca,'FontSize',20);
        end
        
 %   end
    
end

subaxis(1,3,3,'SpacingHoriz',0.075);
Bp=sum(B(:,1:numStrains),2);
Bf=sum(B(:,numStrains+1:end),2);
semilogy(time, Bp,'-', 'Color','k'); hold on;
semilogy(time, Bf,'--', 'Color', 'k'); hold on;
ylabel('Density','fontsize',24);
axis([0 time(end) 1 2*maxY]);
xlabel('Time','fontsize',24);
title(['\gamma=',num2str(params.conj_rate,'%1.0e')]);
if T==24
    xticks(0:T/6:T);
end
yticks(10.^linspace(0,9,10));
set(gca,'FontSize',20);
hleg3=legend({'Plasmid-bearing','Plasmid-free'},'Location','SouthEast');
set(hleg3,'FontSize',20);



%{
if numStrains<6 && tot(end)>0
    
    subaxis(1,3,1);
    hleg1=legend(p_prop,labels,'Location','SouthWest');
    set(hleg1,'FontSize',16);
end
%}

