% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Supplementary Figure 10 in
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');

runsDir='../../data/runs/suppfigure10-11/';
figsDir='../../figures/';

%% CREATE DIRECTORY STRUCTURE

if ~isempty(figsDir)
    if ~exist(figsDir,'dir')
        mkdir(figsDir);
        disp(['mkdir ',figsDir]);
    end
end

%% LOAD DATA

disp([newline,'====== Supp Figure 10-11 ======',newline]);

load([runsDir,'experiment.mat']);
load([runsDir,'sigma_params.mat']);


%%  PLOT SIMULATIONS

simsPath=[runsDir,'sims/'];

which_sims=[1 1 1 1];   %Selected simulation
which_Ms=[1 15 1 15];
which_Cs=fliplr([0 10.^(-10)]);
which_sigmas= [0 0 0.08 0.08];


for isims=1:length(which_sims)
    this_sims=which_sims(isims);
    
    this_sigma=which_sigmas(isims);
    
    this_Ms=which_Ms(isims);
    
    
    plotSimulations(simsPath, [this_sims], [this_sigma], [this_Ms], which_Cs);
    
    rows=length(which_Cs);
    cols=1;
    desc=['sim=',sprintf('%0*d',3,this_sims),' sigma=',num2str(round(this_sigma(1)*100)),'e-2 M=',num2str(which_Ms(isims)),newline,' conjRate=(',num2str(which_Cs,'%1.0e '),')',newline];
    subaxis(rows,cols,1);
    title(desc,'FontSize',20);
    
    if isims==1
        eval(['export_fig ',figsDir,'SuppFigure10de.pdf']);
    elseif isims==2
        eval(['export_fig ',figsDir,'SuppFigure10fg.pdf']);
    elseif isims==3
        eval(['export_fig ',figsDir,'SuppFigure11de.pdf']);
    elseif isims==4
        eval(['export_fig ',figsDir,'SuppFigure11fg.pdf']);
    end
    close();
    
end
