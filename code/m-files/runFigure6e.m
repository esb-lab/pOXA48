% This is the repository for the Matlab codes of the numerical simulations of plasmid dynamics in complex communities.
% The scripts can be used to generate Figure 6e of the article:
% "The distribution of plasmid fitness effects explains plasmid persistence in bacterial communities"
%
% July 27th, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% DEFINE PARAMETERS

maxN=1e4;   % Individuals in meta-community
numExperiments=1e2; %5e3 % Number of simulations

%Checkerboard experiment
Ms=1:1:20;
Cs=1.5e-11;
sigmas=sqrt([0 0.0025 0.01 0.0225 0.004]);
maxSims=100;  % T_sim=maxSims x T
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

cmap_freq=makeColorMap(color_TC, [.95 .95 .95],  color_WT, numColors);

%File structure
mcmcPath='../../data/MCMC_params.mat';
runDir='../../data/runs/figure6e/';

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

%Save experiment parameters
save([runDir,'experiment.mat']);


%% LOAD DATA


disp([newline,'====== run Figure 6e ======',newline]);

%MCMC data
load(mcmcPath);

%Fit bivariate Normal distribution to parameters
gm=getGMModel(MCMC_muKs, MCMC_rhos, max_muKs, max_rhos);

%For data filtering
iTC=find(strcmp(MCMC_plasmids,'TC' ));
iWT=find(strcmp(MCMC_plasmids,'WT' ));

iK=find(strcmp(MCMC_species,'K' ));
iE=find(strcmp(MCMC_species,'E' ));

mean_muKs=[mean(MCMC_muKs(iTC)) mean(MCMC_muKs(iWT))];
mean_rhos=[mean(MCMC_rhos(iTC)) mean(MCMC_rhos(iWT))];
std_muKs=[std(MCMC_muKs(iTC)) std(MCMC_muKs(iWT))];
std_rhos=[std(MCMC_rhos(iTC)) std(MCMC_rhos(iWT))];


%% CREATE COMMUNITIES

sigma_params={};
for si=1:length(sigmas)
    
    %From data:
    cstar=mean_rhos(2);  %Plasmid-free
    Vstar=mean_muKs(2);
    Kstar=1;
    sigma_growth=[std_rhos(2)/max_rhos, std_muKs(2)/max_muKs, 0];  %Varability in c, V, K
    sigma_cost=sigmas(si);  %Variability in cost
    meanCost=.985; %From MCMC parametrization
    
    this_sigma_params=randParams(maxN, T, d, S0, meanCost, cstar, Vstar, Kstar, seg_rate, 0, epsilon, sigma_cost, sigma_growth);
    sigma_params{si}=this_sigma_params;
    
end
save([runDir,'sigma_params.mat'],'sigmas','sigma_params');

%% SIMULATE INVASION EXPERIMENT


mean_pfs=zeros(length(Ms), length(Cs), length(sigmas));
mean_relFitness=zeros(length(Ms), length(Cs), length(sigmas));
for si=1:length(sigmas)
    
    this_sigma=sigmas(si);
    
    disp(['_____ sigma=',num2str(this_sigma)]);
    
    this_std_muKs=std_muKs.*[this_sigma, 1];
    this_std_rhos=std_rhos.*[this_sigma, 1];
    
    % RUN MANY SIMULATIONS OF SUBSETS OF THE COMMUNITY
    all_dyn=zeros(length(Ms), length(Cs), numExperiments);
    all_t_end=zeros(length(Ms), length(Cs), numExperiments);
    all_pf=zeros(length(Ms), length(Cs), numExperiments);
    all_relFitness=zeros(length(Ms), length(Cs), numExperiments);
    
    %par
    for n=1:numExperiments
        
        if mod(n,1e1)==0
            fprintf('.');
            if mod(n,1e2)==0
                disp([num2str(n),' ']);
            end
        end
        
        expe_dyn=zeros(length(Ms), length(Cs));
        expe_t_end=zeros(length(Ms), length(Cs));
        expe_pf=zeros(length(Ms), length(Cs));
        expe_relFitness=zeros(length(Ms), length(Cs));
        for s=1:length(Ms) %Simulate subset
            
            subStrains=Ms(s);
            expe_params=subsetParameters(sigma_params{si}, subStrains);
            expe_params.d=d;
            expe_params.T=T;
            expe_params.colors=cmap_strains(expe_params.index,:);
            
            %ic=getIniConditions(expe_params); %vector of initial conditions
            
            
            B0_WT=zeros(1,subStrains);
            B0_TC=(x0/subStrains).*ones(1,subStrains);
            ic=[S0, B0_TC, B0_WT];
            
            for c=1:length(Cs) %Simulate multiple conjugation rates
                expe_params.conj_rate=Cs(c);
                
                [times, ys, dyn, t_end, pf]=simulateExtinctionMany(expe_params, ic,maxSims); % solve ode
                
                rf=log(sum(ys(end,2:Ms(s))))/log(sum(ys(end,Ms(s)+1:end))); %pair-wise comparison
                
                expe_dyn(s,c)=dyn;
                expe_t_end(s, c)=t_end(1);
                expe_pf(s, c)=pf;
                expe_relFitness(s, c)=rf;
            end
        end
        all_dyn(:,:,n)=expe_dyn;
        all_t_end(:,:,n)=expe_t_end;
        all_pf(:,:,n)=expe_pf;
        all_relFitness(:,:,n)=expe_relFitness;
    end
    %Save file
    save([runDir,'OV_data_sigma',num2str(round(this_sigma*100)),'e-2.mat'],'all_dyn','all_t_end','all_pf','all_relFitness','Ms','Cs','T','d','numExperiments');
    
    this_mean_pf=mean(all_pf(:,:,:), 3);
    mean_pfs(:,:,si)=this_mean_pf;
    
    this_mean_rf=mean(all_relFitness(:,:,:),3);
    mean_relFitness(:,:,si)=this_mean_rf;
    
end
