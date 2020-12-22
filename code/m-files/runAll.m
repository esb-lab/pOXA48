% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate the figures in the main text and supplementary information of the article
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clc
close all
clear all
warning off export_fig:exportgraphics

%% FILE STRUCTURE

figDir='../../figures/';
if ~exist(figDir,'dir')
    mkdir(figDir);
    disp(['mkdir ',figDir]);
end

%%


%Figure 5a-b: Modelling pOXA-48 fitness effects.
plotFigure5ab;

%Figure 5c: Fraction of plasmid-bearing as a function of conjugation rate
%runFigure5c;
plotFigure5c;

%Figure 5d: Fraction of plasmid-bearing cells as a function of the rate of HGT
%runFigure5d;
plotFigure5d;

%Figure 6: Modelling plasmid persistence in polymicrobial communities
%runFigure6;
plotFigure6;

%Figure 6e: Mean fraction of plasmid-bearing cells as a function of the number of strains in the community
%runFigure6e;
plotFigure6e;

%Supp Figure 7:In silico competition experiments
plotSuppFigure7

%Supp Figure 8: Comparison between in silico and experimental competition experiments
plotSuppFigure8

%Supp Figure 9. Plasmid stability as a function of the conjugation rate and plasmid cost
plotSuppFigure9;

%Supp Figure 10: Effect of conjugation and community complexity in plasmid dynamics
%runSuppFigure10;
plotSuppFigure10;

%Supp Figure 15: Sampling the distribution of plasmid fitness effects
plotSuppFigure15;

%Supp File 1. MCMC Diagnostic plots
plotSuppFile1;
runGelmanRubin;