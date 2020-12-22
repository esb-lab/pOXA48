% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Supplementary Figure 9 in
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% PARAMETERS


disp([newline,'====== Supp Figure 9 ======',newline]);

expeDir='../../data/runs/figure5d/';
figDir='../../figures/';


load([expeDir,'experiment.mat']);  %Parameters

%%
%
this_sigma=0.000;
figure();
clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);

p=[];

sel_mus=linspace(0.0, 0.5, 11);
%sel_mus=mus;
cmap_mu=cbrewer('seq', 'Blues', length(sel_mus)+5);
cmap_mu=cmap_mu(5:end,:);
for mi=1:length(sel_mus)
    
    
    matFile=[expeDir,'conjRate_assay_mu',num2str(round(sel_mus(mi)*1e3)),'e-3_sigma',num2str(round(this_sigma*1e3)),'e-3.mat'];
    if exist(matFile, 'file')
        disp(['Loading ',matFile]);
        load(matFile);
        
        this_style='-';
        this_color=cmap_mu(mi,:);
        
        semilogx(conj_rates, 50.*ones(1,length(conj_rates)),'k:'); hold on;
        p(mi)=semilogx(conj_rates, 100*smooth(expe_pFraction(:,1)), this_style,'LineWidth',2, 'Color',this_color);
    else
        disp(['Not found ',matFile]);
    end
end
%
colormap(cmap_mu);
cbh=colorbar;
yt=linspace(0,1,length(sel_mus)+1);
set(cbh,'YTick',yt(1:end-1)+yt(2)/2)
set(cbh,'YTickLabel',100*(sel_mus))
title(cbh, 'Cost (%)')
%
xlim([0 max(conj_rates)]);
ylim([0 100]);
set(gca,'FontSize',20);
xlabel('Conjugation rate (\gamma)','FontSize',24);
ylabel('Plasmid fraction (%)','FontSize',24);

eval(['export_fig ',figDir,'SuppFigure9.pdf']);

