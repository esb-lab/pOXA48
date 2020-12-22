% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Figure 6 in 
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% CREATE DIRECTORY STRUCTURE

expeDir='../../data/runs/figure6/';
outDir='../../figures/';

%% LOAD DATA

load([expeDir,'experiment.mat']);  %Parameters
load([expeDir,'sigma_params.mat']);

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


%% LOAD COMMUMITIES

dataFile='sigma_params.mat';
pathFile=[expeDir,dataFile];
if exist(pathFile,'file')
    load(pathFile);
    disp(['Loading ',pathFile]);
else
    disp(['Not found ',pathFile]);
end

%% PLOT COMMUNITIES

for si=1:length(sigma_params)
    figure(); clf('reset');
    set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white');
    set(gcf,'units','points','position',[616   598   560/2   420/2])
    
    plotGMModel(sigma_params{si}, max_muKs, max_rhos);
    
    if sigmas(si)==0
        eval(['export_fig ',outDir,'Figure6a_inset.pdf']);
    elseif sigmas(si)==0.0837  %Data
        eval(['export_fig ',outDir,'Figure6b_inset.pdf']);
    else
        eval(['export_fig ',outDir,'Figure6_inset_sigma',num2str(sigmas(si)*100),'e-2.pdf']);
    end
end

%% Load DFE

dataFile='DFE_sigmas.mat';
pathFile=[expeDir,dataFile];
if exist(pathFile,'file')
    load(pathFile);
    disp(['Loading ',pathFile]);
else
    disp(['Not found ',pathFile]);
end

%% plotDFE

disp([newline,'====== Figure 6ab ======',newline]);

if exist('sigma_relFitness', 'var')
    
    nbins=41;
    bin_min=.75;
    bin_max=1.25;
    bin_delta_sample=.025;
    
    for is=1:length(sigmas)
        
        figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);
        this_relFitness=sigma_relFitness(is,:);
        
        ebins_sample=(bin_min:bin_delta_sample:bin_max)+-(1-mean(this_relFitness))/2;
        
        h_sample=histogram(this_relFitness,ebins_sample,  'FaceColor',[.9 .9 .9], 'FaceAlpha',1,'EdgeColor','w', 'Normalization','count');
        sample_y=h_sample.Values;
        
        sample_x=ebins_sample(2:end)-bin_delta_sample/2;
        
        disp(['sigma=',num2str(sigmas(is)),newline,'   Mean(gm_relFitness)=',num2str(mean(this_relFitness))]);
        disp(['   Std(gm_relFitness)=',num2str(std(this_relFitness))]);
        
        bar(sample_x, sample_y ,'LineWidth',1,'FaceColor',[.9 .9 .9]); hold on;
        plot([mean(this_relFitness) mean(this_relFitness)],[0 1.1*max(sample_y)],'--r','LineWidth',2);
        
        plot([1 1],[0 1.1*max(sample_y)],'-k','LineWidth',1);
        xlim([0.7 1.5]);
        ylim([0 1.1*max(sample_y)]);
        xticks(0:.2:2);
        
        set(gca,'FontSize',20);
        xlabel('Relative fitness','FontSize',24);
        ylabel('Count','FontSize',24);
        
        title(['\sigma=',num2str(sigmas(is))]);
        
        if sigmas(is)==0
            eval(['export_fig ',outDir,'Figure6a.pdf']);
        elseif sigmas(is)==0.0837  %Data
            eval(['export_fig ',outDir,'Figure6b.pdf']);
        else
            eval(['export_fig ',outDir,'Figure6_sigma',num2str(sigmas(si)),'.pdf']);
        end
    end
else
    disp('DFE data not found');
end


%% PLOT INVASION EXPERIMENT

disp([newline,'====== Figure 6cd ======',newline]);

mean_pfs=zeros(length(Ms), length(Cs), length(sigmas));
mean_relFitness=zeros(length(Ms), length(Cs), length(sigmas));
std_pfs=zeros(length(Ms), length(Cs), length(sigmas));
std_relFitness=zeros(length(Ms), length(Cs), length(sigmas));
for si=1:length(sigmas)
    this_sigma=sigmas(si);
    
    dataFile=['OV_data_sigma',num2str(round(this_sigma*100)),'e-2.mat'];
    pathFile=[expeDir,dataFile];
    if exist(pathFile,'file')
        
        %Load data file
        load(pathFile);
        disp(['Loading ',pathFile]);
        
        %Plot simulations
        
        plotPlasmidFractionCheckerboard;
        if sigmas(si)==0
            eval(['export_fig ',outDir,'Figure6c.pdf']);
        elseif sigmas(si)==0.0837  %Data
            eval(['export_fig ',outDir,'Figure6d.pdf']);
        else
            eval(['export_fig ',outDir,'Figure6_pf_sigma',num2str(round(this_sigma*100)),'e-2.png']);
        end
        
        %Store statistics
        this_mean_pf=mean(all_pf(:,:,:), 3);
        mean_pfs(:,:,si)=this_mean_pf;
        
        this_mean_rf=mean(all_relFitness(:,:,:),3);
        mean_relFitness(:,:,si)=this_mean_rf;
        
        this_std_pf=std(all_pf(:,:,:), 0, 3);
        std_pfs(:,:,si)=this_std_pf;
        
        this_std_rf=std(all_relFitness(:,:,:),0,3);
        std_relFitness(:,:,si)=this_std_rf;
        
    else
        disp(['Not found: ',pathFile]);
    end
end

