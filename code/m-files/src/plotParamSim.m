function plotParamSim(time, y, params, max_muKs, max_rhos, gm, dyn)


if nargin<3
    strain_labels={};
    strain_colors=jet(params.numStrains);
else
    if isfield(params,'strains')
        strain_labels=params.code;
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
level_color=[.9 .9 .9];
rmax=2400;
rmin=24;


%figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'fontsize',20);
%set(gcf, 'Units','normalized','Position',[0. 0.2 0.6 .3]);

%T=time(end);

%S=y(:,1);
B=y(:,2:end);


[~, totStrains]=size(B);
numStrains=totStrains/2;

tot=sum(B,2);  %Total density of bacteria

maxY=max([max(sum(B(:,1:numStrains),2)), max(sum(B(:,numStrains+1:end),2)) ]);
labels={};
p_prop=zeros(1,numStrains);

% Plot contour of metacommunity
%subaxis(1,3,1,'SpacingHoriz',0.075,'PaddingBottom',.1);
if ~isempty(gm)
    fcontour(@(x,y)reshape(pdf(gm,[x(:),y(:)]),size(x)),[0 1],'LineWidth',.5,'LineColor',level_color,'MeshDensity',100); hold on;
end


for i=1:numStrains
    
%    if ~strcmp(strain_labels{i},'X')

        
        if tot(end)>0

            this_color=strain_colors(i,:);

            %Plot parameters w/ circle radius proportional to relative density
            %subaxis(1,3,1,'SpacingHoriz',0.075,'PaddingBottom',.1);
            
            this_freqTC=B(end,i)./tot(end);
            this_freqWT=B(end,i+params.numStrains)./tot(end);

            rmarker_TC=rmin+(rmax-rmin)*this_freqTC;
            rmarker_WT=rmin+(rmax-rmin)*this_freqWT;

            %Plasmid-bearing
            %scatter(params.Vs(i)/max_muKs, params.cs(i)/max_rhos,rmarker_WT,this_color,'filled'); hold on;
            s_TC=scatter(params.Vs(i)/max_muKs, params.cs(i)/max_rhos,rmarker_WT,this_color); hold on;
            s_TC.LineWidth = 1;
            s_TC.MarkerEdgeColor = this_color;
            s_TC.MarkerFaceColor = this_color;
            s_TC.MarkerFaceAlpha = 0.5;
            
            %Plasmid-free
            s_WT=scatter(params.Vs(i+numStrains)/max_muKs, params.cs(i+numStrains)/max_rhos, rmarker_TC,this_color,'MarkerFaceAlpha',1); hold on;
            s_WT.LineWidth = 2;
            s_WT.MarkerEdgeColor = this_color;
            
            %plot(params.Vs(i+params.numStrains)/max_muKs,params.cs(i+params.numStrains)/max_rhos, str_marker,'LineWidth',2,'MarkerEdgeColor',this_color,'MarkerSize',rmarker_WT); hold on;
            %labels{i}=extractBefore(strain_labels{i},'_');

            set(gca,'fontsize',20);
            %title(['M=',num2str(params.numStrains)],'FontSize',20);
            xlabel('V_{max}/K_m x 10^{-10}','FontSize',24);
            ylabel('\rho x 10^{8}','FontSize',24);
            nticks=5;
            axis([0. 1 0. 1]);
            xticks(0:1/nticks:1);
            xticklabels((0:max_muKs/nticks:max_muKs)/1e-9);
            yticks(0:1/nticks:1);
            yticklabels((0:max_rhos/nticks:max_rhos)/1e9);
            text(0.1,0.1,strain_labels{i},'FontSize',24,'Color',this_color)
            
        end
        
 %   end
    
end

