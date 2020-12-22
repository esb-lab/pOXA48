% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate the figures in the main text and supplementary information of the article
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx

clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');

color_TC=[0     0.56863     0.56863];
color_WT=[0.83137     0.36863           0];

%% PARAMETERS

% Open output files from MCMC and load posterior distributions
rootDir = '../../data/MCMC_chains_gr/';

all_codes={ 'Ec01',  'Ec02',  'Ec03',  'Ec04',  'Ec05',  'Ec06',  'Ec07',  'Ec08',  'Ec09',  'Ec10',  'Ec11',  'Ec12',  'Ec13',  'Ec14',  'Ec15',  'Ec16',  'Ec17',  'Ec18',  'Ec19',  'Ec20',  'Ec21',  'Ec22',  'Ec23',  'Ec24',  'Ec25',  'Kpn01',  'Kpn02',  'Kpn03',  'Kpn04',  'Kpn05',  'Kpn06',  'Kpn07',  'Kpn08',  'Kpn09',  'Kpn10',  'Kpn11',  'Kpn13',  'Kpn14',  'Kpn15',  'Kpn16',  'Kpn17',  'Kpn18',  'Kpn19',  'Kpn20',  'Kpn21',  'Kpn22',  'Kpn23',  'Kpn24',  'Kpn25'};
all_strains={ 'C001',  'C002',  'C006',  'C011',  'C012',  'C021',  'C022',  'C031',  'C051',  'C063',  'C094',  'C107',  'C115',  'C131',  'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',  'K141',  'K168',  'K177',  'K200',  'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
all_plasmids={'WT', 'TC'};

M=3; %Number of chains

%% Compute Gelman-Rubin diagnostic for all strains

for istrain=1:length(all_strains)
    for itype=1:length(all_plasmids)

        lbl_strain=all_strains{istrain};
        lbl_code=all_codes{istrain};
        lbl_prior='uniform';
        lbl_clones='nclones1';
        lbl_type=all_plasmids{itype};

        dataDir=[rootDir,lbl_strain,'_',lbl_prior,'_',lbl_clones,'_',lbl_type,'/'];
        
        [muK_R, rho_R]=computeGelmanRubin(dataDir, M);
        
        disp([lbl_code,'_',lbl_type,': R_muK=',num2str(muK_R),' | R_rho=',num2str(rho_R)])
        
    end
end