% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Supp Figure 15 in the article
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% PARAMETERS

T=24;   %Length of experiment
S0=1;    %Resource concentration
x0=1e6;  %Initial bacterial density
d=0.0;  %dilution rate
epsilon=1e-6; %extinction threshold
seg_rate=1e-8;  %Segregation rate

max_muKs=10e-10; %Limits for gmmodel
max_rhos=12e8;

color_TC=[0     0.56863     0.56863];
color_WT=[0.83137     0.36863           0];

mcmcPath='../../data/MCMC_params.mat';
figDir='../../figures/';

%% SIMULATION PARAMETERS

%From data
strains={'C001', 'C022', 'C002',  'C006',  'C011',  'C012',  'C021',  'C031',  'C051',  'C063',    'C107',  'C115',   'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',    'K168',  'K177', 'K200',    'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
plasmids={'WT','TC'};
totStrains=length(strains);  % Number of strains (data)

%For conjugation assay
numExperiments=10;
subStrains=2; %Size of sub-community

conj_rates=[1e-9];%logspace(-11, -8, 11); %Conjugation rates used to determine optimal

mus=linspace(0.0, 0.25, 5);
sigmas=sqrt([0 0.0075 0.015 0.0225]);  %DFE standard deviation

%% PATHS AND FILE STRUCTURE

runDir='runs/suppfigure15/';
if ~exist(runDir, 'dir')
    mkdir(runDir);
    disp(['mkdir ',runDir]);
end

%Colormap
cmap_brewer=cbrewer('qual', 'Paired', totStrains*2);
cmap=cmap_brewer([1:2:totStrains*2 2:2:totStrains*2], :);


%% LOAD DATA

%MCMC data
load(mcmcPath);

%For data filtering
iTC=find(strcmp(MCMC_plasmids,'TC' ));
iWT=find(strcmp(MCMC_plasmids,'WT' ));

iK=find(strcmp(MCMC_species,'K' ));
iE=find(strcmp(MCMC_species,'E' ));

mean_muKs=[mean(MCMC_muKs(iTC)) mean(MCMC_muKs(iWT))];
mean_rhos=[mean(MCMC_rhos(iTC)) mean(MCMC_rhos(iWT))];

std_muKs=[std(MCMC_muKs(iTC)) std(MCMC_muKs(iWT))];
std_rhos=[std(MCMC_rhos(iTC)) std(MCMC_rhos(iWT))];

%% Find critical conjugation rate for different mus

disp([newline,'====== Supp Figure 15 ======']);

figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white');
set(gcf, 'Units','normalized','Position',[0. 0. 1 1]);

for mi=1:length(mus)
    
    mean_cost=1-mus(mi);
    
    for si =1:length(sigmas)
        
        %Variability in cost
        sigma_cost=sigmas(si);
        
        for icr=1:length(conj_rates)
            
            conj_rate=conj_rates(icr);
            
            relFitness=zeros(1,numExperiments);
            pFractions=zeros(1,numExperiments);
            for n=1:numExperiments
                
                %Parameters
                this_colors=[rand rand rand];
                cstar=mean_rhos(2);
                Vstar=mean_muKs(2);
                Kstar=1;
                sigma_growth=[std_rhos(2)/max_rhos, std_muKs(2)/max_muKs, 0];  %Varability in c, V, K
                
                this_params=randParams(subStrains, T, d, S0, mean_cost, cstar, Vstar, Kstar, seg_rate, conj_rate, epsilon, sigma_cost, sigma_growth);
                
                %Plot community
                subaxis(length(sigmas), length(mus),mi,length(sigmas)- si+1, 'SpacingVert',0.05,'SpacingHoriz',0.05);
                plotParams(this_params, max_muKs, max_rhos, si)
                title(['w=',num2str(1-mus(mi)),', \sigma^2=',num2str(sigma_cost^2)]);
                if si==1 && mi==1
                    
                elseif si==1 && mi>1
                    set(gca,'YTickLabel','');
                    ylabel('');
                elseif mi==1
                    set(gca,'XTickLabel','');
                    xlabel('');
                else
                    xlabel('');
                    ylabel('');
                    set(gca,'YTickLabel','');
                    set(gca,'XTickLabel','');
                end
            end
        end
    end
end

eval(['export_fig ',figDir,'SuppFigure15.pdf']);
close();
