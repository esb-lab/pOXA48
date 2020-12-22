% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Supplementary File 1 (MCMC Diagnostic tests)
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


%% PARAMETERS
all_codes={ 'Ec01',  'Ec02',  'Ec03',  'Ec04',  'Ec05',  'Ec06',  'Ec07',  'Ec08',  'Ec09',  'Ec10',  'Ec11',  'Ec12',  'Ec13',  'Ec14',  'Ec15',  'Ec16',  'Ec17',  'Ec18',  'Ec19',  'Ec20',  'Ec21',  'Ec22',  'Ec23',  'Ec24',  'Ec25',  'Kpn01',  'Kpn02',  'Kpn03',  'Kpn04',  'Kpn05',  'Kpn06',  'Kpn07',  'Kpn08',  'Kpn09',  'Kpn10',  'Kpn11',  'Kpn13',  'Kpn14',  'Kpn15',  'Kpn16',  'Kpn17',  'Kpn18',  'Kpn19',  'Kpn20',  'Kpn21',  'Kpn22',  'Kpn23',  'Kpn24',  'Kpn25'};
all_strains={ 'C001',  'C002',  'C006',  'C011',  'C012',  'C021',  'C022',  'C031',  'C051',  'C063',  'C094',  'C107',  'C115',  'C131',  'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',  'K141',  'K168',  'K177',  'K200',  'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
all_plasmids={'WT', 'TC'};

T=24;   %Length of experiment
S0=1;    %Resource concentration

d=0.0;  %dilution rate

seg_rate=1e-8;  %Segregation rate
conj_rate=0e-12;  %Conjugation rate

epsilon=1e-6; %extinction threshold

max_muKs=10e-10; %Limits for gmmodel
max_rhos=12e8;

color_WT=[0     0.56863     0.56863];
color_TC=[0.83137     0.36863           0];

%Colormap
cmap_brewer=cbrewer('qual', 'Paired', length(all_strains)*2);
cmap=cmap_brewer([1:2:length(all_strains)*2 2:2:length(all_strains)*2], :);

%% PATHS

figDir='../../figures/SuppFile1/';
rootDir = '../../data/MCMC_chains/';
strainPath='../../data/ODdata_strains/';

disp([newline,'====== Supp File 1 ======',newline]);

%% LOAD EXPERIMENTAL DATA

%Flow cytometry data (for comparison)
DATA_strains={'C001', 'C002', 'C006', 'C011', 'C012', 'C021', 'C022', 'C031', 'C051', 'C063', 'C094', 'C107', 'C115', 'C131', 'C141', 'C201', 'C227', 'C232', 'C247', 'C261', 'C286', 'C290', 'C302', 'C309', 'C324', 'K037', 'K038', 'K087', 'K094', 'K112', 'K114', 'K125', 'K141', 'K168', 'K177', 'K200', 'K201', 'K209', 'K213', 'K216', 'K224', 'K225', 'K241', 'K248', 'K249', 'K253', 'K257', 'K275', 'K285', 'K300'};
DATA_relFitness=[0.924055921	1.005585204	0.953955342	0.8697700094	0.9790679309	0.8612088774	0.9836714827	0.976027132	1.005880695	1.025256181	1.03845763	1.05893746	1.001274491	1.151862768	1.0642684	1.03421127	1.00242249	0.8617113414	0.9804258966	0.9873267401	0.8784454426	1.01482229	0.9905312712	0.8253910816	0.8103264669	0.9170508336	0.985061177	0.9928021947	0.9141655888	0.9059386546	0.9560160898	1.014466675	0.7667828255	1.015722835	1.036012759	0.8967216527	0.9658668826	0.7616641848	0.9445656777	1.17317741	0.890155927	0.9761680929	1.023776434	1.043786826	0.9496171203	1.051516039	1.134035406	0.9715462481	0.9826900563	0.9785162211];

mcmcPath='../../data/MCMC_params.mat';

%MCMC data
load(mcmcPath);

%% CREATE DIRECTORY STRUCTURE

if ~isempty(figDir)
    if ~exist(figDir,'dir')
        mkdir(figDir);
        disp(['mkdir ',figDir]);
    end
end

selected_strains=all_strains;
%selected_strains={'C001' };  %Comment this line to produce all plots

%%
for s=1:length(selected_strains)
    
    
    figure(); clf('reset');
    set(gcf, 'Units','normalized','Position',[0. 0. .6 1]);
    set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); hold on;
    
    %PRIORS
    plotMCMCpriors;
    
    %REALIZATIONS
    subaxis(4, 4, 3, 2,2,2,'PaddingTop',0.05,'PaddingBottom',0.05,'SpacingHoriz',0.0);
    plotMCMCdensity;
    
    %DATA CLONING
    plotMCMCclones;
    
    % CHAINS
    plotMCMCchains;
    plotMCMCchains;
    
    %POSTERIOR
    subaxis(4, 4, 1, 3, 2, 2,'PaddingTop',0.05,'PaddingRight',0.05);
    plotMCMCposterior;
    
    eval(['export_fig ',figDir,'MCMC_',code_strain,'.pdf']);
    close;
    
    disp(['Exporting ',selected_strains{s},': ','MCMC_',code_strain,'.pdf']);
    
end
