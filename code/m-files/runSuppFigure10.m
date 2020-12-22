% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Supp Figure 10 and Supp Figure 11 of the  article
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% DEFINE PARAMETERS

% Individuals in meta-community
maxN=1e3;  
numExperiments=10; %
startExperiment=1;

%Checkerboard experiment
Ms=[1 15]; 
Cs=[0 10.^(-12:2:-8)];
sigmas= [0, 0.0837];

maxSims=10;  % T_sim=maxSims x T
T=100;   %Maximum length of experiment (interval)
S0=1;    %Resource concentration
x0=1e6;  %Initial bacterial density

d=0.2;  %dilution rate

seg_rate=1e-8;  %Segregation rate
epsilon=1; %extinction threshold  

max_muKs=10e-10; %Limits for gmmodel
max_rhos=12e8;

%Color scheme
color_WT=[0     0.56863     0.56863];
color_TC=[0.83137     0.36863           0];
numColors=201;

cmap = distinguishable_colors(numColors);
icolors=repmat(1:numColors, 1, round(maxN/numColors));
icolors=icolors(randperm(length(icolors)));
icolors=icolors(1:maxN);

cmap_strains=zeros(maxN, 3);
for i=1:maxN
   cmap_strains(i,:)=cmap(icolors(i),:);
end

cmap_freq=makeColorMap(color_WT, [.95 .95 .95],  color_TC, numColors);

%File structure
mcmcPath='../../data/MCMC_params.mat';
runDir='../../data/runs/suppfigure10-11/';

simDir=[runDir,'sims/'];  %Directory to save realizations


%% CREATE DIRECTORY STRUCTURE

if ~isempty(runDir)
    if ~exist(runDir,'dir')
        mkdir(runDir);
        disp(['mkdir ',runDir]);
    end
end

if ~exist(runDir,'dir')
        mkdir(runDir);
        disp(['mkdir ',runDir]);
end

if ~isempty(simDir)
    if ~exist(simDir,'dir')%
        mkdir(simDir);
        disp(['mkdir ',simDir]);
        mkdir([simDir,'PDF/']);
        disp(['mkdir ',simDir,'PDF/']);
        mkdir([simDir,'PNG/']);
        disp(['mkdir ',simDir,'PNG/']);
    end
end

%Save experiment parameters
save([runDir,'experiment.mat']);


%% LOAD DATA

disp([newline,'====== LOADING DATA ======',newline]);

%MCMC data
load(mcmcPath);
disp([num2str(length(MCMC_strains)),' MCMC output files loaded']);

%Fit bivariate Normal distribution to parameters
gm=getGMModel(MCMC_muKs, MCMC_rhos, max_muKs, max_rhos);

%For data filtering
iTC=find(strcmp(MCMC_plasmids,'TC' ));
iWT=find(strcmp(MCMC_plasmids,'WT' ));

iK=find(strcmp(MCMC_species,'K' ));
iE=find(strcmp(MCMC_species,'E' ));

%Fit bivariate Normal distribution to TC/WT parameters
gm_TC=getGMModel(MCMC_muKs(iTC), MCMC_rhos(iTC), max_muKs, max_rhos);
gm_WT=getGMModel(MCMC_muKs(iWT), MCMC_rhos(iWT), max_muKs, max_rhos);

mean_muKs=[mean(MCMC_muKs(iTC)) mean(MCMC_muKs(iWT))];
mean_rhos=[mean(MCMC_rhos(iTC)) mean(MCMC_rhos(iWT))];

std_muKs=[std(MCMC_muKs(iTC)) std(MCMC_muKs(iWT))];
std_rhos=[std(MCMC_rhos(iTC)) std(MCMC_rhos(iWT))];

data_gm_TC = gmdistribution([mean_muKs(1) mean_rhos(1)],[std_muKs(1) std_rhos(1)]);
data_gm_WT = gmdistribution([mean_muKs(2) mean_rhos(2)],[std_muKs(2) std_rhos(2)]);

%% CREATE COMMUNITIES

sigma_params={};
for si=1:length(sigmas)
    
    %From data:
    cstar=mean_rhos(2);
    Vstar=mean_muKs(2); 
    Kstar=1;
    sigma_growth=[std_rhos(2)/max_rhos, std_muKs(2)/max_muKs, 0];  %Varability in c, V, K
    sigma_cost=sigmas(si);  %Variability in cost
    meanCost=.95;
    
    this_sigma_params=randParams(maxN, T, d, S0, meanCost, cstar, Vstar, Kstar, seg_rate, 0, epsilon, sigma_cost, sigma_growth);
    sigma_params{si}=this_sigma_params;
    
end
save([runDir,'sigma_params.mat'],'sigmas','sigma_params');


%% SIMULATE INVASION EXPERIMENT (and plot state variables)


disp([newline,'====== SIMULATE MULTI-STRAIN EXPERIMENT ======',newline]);
for si=1:length(sigmas)

    this_sigma=sigmas(si);
    
    disp(['_____ sigma=',num2str(this_sigma)]);
     
    this_std_muKs=std_muKs.*[this_sigma, 1];
    this_std_rhos=std_rhos.*[this_sigma, 1];
     
    % RUN MANY SIMULATIONS OF SUBSETS OF THE COMMUNITY
    for n=startExperiment:numExperiments
        
        tic
        for s=1:length(Ms) %Simulate subset
            
            subStrains=Ms(s);
            expe_params=subsetParameters(sigma_params{si}, subStrains);
            expe_params.d=d;
            expe_params.T=T;
            expe_params.colors=cmap_strains(expe_params.index,:);
            
            %Initial conditions
            B0_WT=zeros(1,subStrains);
            B0_TC=(x0/subStrains).*ones(1,subStrains);
            ic=[S0, B0_TC, B0_WT];
            
            for c=1:length(Cs) %Simulate multiple conjugation rates
                expe_params.conj_rate=Cs(c);

                [times, ys, dyn, t_end, pf]=simulateExtinctionMany(expe_params, ic,maxSims); % solve ode
                
                %Plot realizations
                if ~isempty(simDir)
                    plotManyGrowthCurves(times, ys, expe_params, max_muKs, max_rhos, gm, dyn);
                    eval(['export_fig ',simDir,'PNG/sim',sprintf('%0*d',3,n),'_sigma',num2str(round(this_sigma*100)),'e-2_M',num2str(Ms(s)),'_conjRate',num2str(Cs(c),'%1.0e'),'.png']);
                    eval(['export_fig ',simDir,'PDF/sim',sprintf('%0*d',3,n),'_sigma',num2str(round(this_sigma*100)),'e-2_M',num2str(Ms(s)),'_conjRate',num2str(Cs(c),'%1.0e'),'.pdf']);
                    close();
                end

            end
            toc
        end
    end
    
end



         