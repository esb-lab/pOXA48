% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Figure 5c in the manuscript
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% PARAMETERS

expeDir='../../data/runs/figure5c/';
figDir='../../figures/';


%%


disp([newline,'====== Figure 5c ======',newline]);

figure();
clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);
p=[];
for si=1:2
    
    if si==1
       load([expeDir,'conjRate_assay_sigmaData.mat']) 
       this_style='-k';
    else
       load([expeDir,'conjRate_assay_sigma0.mat'])
       this_style='--k';
    end
    
    semilogx(conj_rates, 50.*ones(1,length(conj_rates)),'k:'); hold on;
    p(si)=semilogx(conj_rates, 100*expe_pFraction(:,1), this_style,'LineWidth',3);
end

xlim([0 max(conj_rates)]);
ylim([0 100]);
set(gca,'FontSize',20);
xlabel('Conjugation rate','FontSize',24);
ylabel('Plasmid fraction (%)','FontSize',24);

eval(['export_fig ',figDir,'Figure5c.pdf']);