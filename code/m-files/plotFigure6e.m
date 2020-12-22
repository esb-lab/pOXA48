% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Figure 6e in the article
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');

expeDir='../../data/runs/figure6e/';
figDir='../../figures/';

%% LOAD DATA

disp([newline,'====== Figure 6e ======',newline]);
load([expeDir,'sigma_params.mat']);
load([expeDir,'experiment.mat']);

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

dataFile=['sigma_params.mat'];
pathFile=[expeDir,dataFile];
if exist(pathFile,'file')
    load(pathFile);
    disp(['Loading ',pathFile]);
else
    disp(['Not found ',pathFile]);
end

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

%% COMPARE DIFFERENT SIGMAS

sigma_colors=makeColorMap(color_WT,color_TC,length(sigmas));
for iCs=1:length(Cs)
    
    leg={};
    
    %Plasmid frequency
    figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);
    set(gcf,'units','points','position',[616   598   560*1.5   420])
    for isc=1:length(sigmas)
        
        this_cs_pf=mean_pfs(:,iCs,isc);
        this_color=sigma_colors(isc,:);
        
        if ~isempty(this_cs_pf)
            plot(Ms(1:end), 100*this_cs_pf(1:end), '-', 'Color', this_color, 'LineWidth',3); hold on;
            
            leg{isc}=['\sigma^2=', num2str(sigmas(isc)^2)];
            
        end
    end
    if ~isempty(leg)
        legend(leg, 'Location','SouthEast', 'FontSize',20);
    end
    set(gca,'FontSize',20)
    xlabel('Number of strains in community', 'FontSize',24);
    ylabel('Plasmid frequency (%)', 'FontSize',24);
    ylim([0 105]);
    xlim([1-0.5 Ms(end)+.5]);
    
    eval(['export_fig ',figDir,'Figure6e.pdf']);
end


