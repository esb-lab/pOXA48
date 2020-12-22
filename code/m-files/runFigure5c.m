% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate the data necessary to produce Figure 5c of the manuscript
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


disp([newline,'====== run Figure 5c ======',newline]);

%% SIMULATION PARAMETERS

strains={'C001', 'C022', 'C002',  'C006',  'C011',  'C012',  'C021',  'C031',  'C051',  'C063',    'C107',  'C115',   'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',    'K168',  'K177', 'K200',    'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
plasmids={'WT','TC'};

totStrains=length(strains);  % Number of strains (data)
simStrains=1e4;  %Number of strains (sample gmmodel)

%For conjugation assay
numExperiments=1e2;  %1e4;
subStrains=2; %Size of sub-community
conj_rates=logspace(-11, -8, 11);  %logspace(-11, -8, 101);

%% PATHS AND FILE STRUCTURE

runDir='../../data/runs/figure5c/';
if ~exist(runDir, 'dir')
    mkdir(runDir);
    disp(['mkdir ',runDir]);
end

mcmcPath='../../data/MCMC_params.mat';

%% LOAD DATA

%Flow cytometry data (for comparison)
DATA_strains={'C001', 'C002', 'C006', 'C011', 'C012', 'C021', 'C022', 'C031', 'C051', 'C063', 'C094', 'C107', 'C115', 'C131', 'C141', 'C201', 'C227', 'C232', 'C247', 'C261', 'C286', 'C290', 'C302', 'C309', 'C324', 'K037', 'K038', 'K087', 'K094', 'K112', 'K114', 'K125', 'K141', 'K168', 'K177', 'K200', 'K201', 'K209', 'K213', 'K216', 'K224', 'K225', 'K241', 'K248', 'K249', 'K253', 'K257', 'K275', 'K285', 'K300'};
DATA_relFitness=[0.924055921	1.005585204	0.953955342	0.8697700094	0.9790679309	0.8612088774	0.9836714827	0.976027132	1.005880695	1.025256181	1.03845763	1.05893746	1.001274491	1.151862768	1.0642684	1.03421127	1.00242249	0.8617113414	0.9804258966	0.9873267401	0.8784454426	1.01482229	0.9905312712	0.8253910816	0.8103264669	0.9170508336	0.985061177	0.9928021947	0.9141655888	0.9059386546	0.9560160898	1.014466675	0.7667828255	1.015722835	1.036012759	0.8967216527	0.9658668826	0.7616641848	0.9445656777	1.17317741	0.890155927	0.9761680929	1.023776434	1.043786826	0.9496171203	1.051516039	1.134035406	0.9715462481	0.9826900563	0.9785162211];

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

mean_muKs=[mean(MCMC_muKs(iTC)) mean(MCMC_muKs(iWT))];
mean_rhos=[mean(MCMC_rhos(iTC)) mean(MCMC_rhos(iWT))];
disp(['Mean muKs (TC,WT)=(',num2str(mean_muKs(1)),',',num2str(mean_muKs(2)),')']);
disp(['Mean rhos (TC,WT)=(',num2str(mean_rhos(1)),',',num2str(mean_rhos(2)),')']);

std_muKs=[std(MCMC_muKs(iTC)) std(MCMC_muKs(iWT))];
std_rhos=[std(MCMC_rhos(iTC)) std(MCMC_rhos(iWT))];
disp(['Std muKs (TC,WT)=(',num2str(std_muKs(1)),',',num2str(std_muKs(2)),')']);
disp(['Std rhos (TC,WT)=(',num2str(std_rhos(1)),',',num2str(std_rhos(2)),')']);

%Colormap
cmap_brewer=cbrewer('qual', 'Paired', totStrains*2);
cmap=cmap_brewer([1:2:totStrains*2 2:2:totStrains*2], :);

%% EVALUATING THE ROLE OF HORIZONTAL TRANSMISSION (SAMPLE)

disp([newline,'====== CONJUGATION ASSAY ======',newline]);

for si=1:2
    
    expe_pFraction=zeros(length(conj_rates),2);
    expe_relFitness=zeros(length(conj_rates),2);

    for icr=1:length(conj_rates)

        disp([newline, num2str(icr),' ** conj_rate=',num2str(conj_rates(icr))]);

        conj_rate=conj_rates(icr);

        relFitness=zeros(1,numExperiments);
        pFractions=zeros(1,numExperiments);
        for n=1:numExperiments

            if mod(n,1e2)==0
                fprintf('.');
                if mod(n,1e3)==0
                    disp([num2str(n),' ']);
                end
            end

            %Random community
            istrains=randperm(length(strains),subStrains); % randi(length(strains), [1, subStrains]);

            this_strains=strains(istrains);
            
            %Parameter values
            this_colors=cmap(istrains,:);
            
            if(si==1)    % SAMPLE RANDOM STRAINS
                this_params=dataParameters(mcmcPath, this_strains, this_colors, T, d, S0, seg_rate, conj_rate, epsilon);
            else
                cstar=mean_rhos(2);
                Vstar=mean_muKs(2); 
                Kstar=1;
                sigma_growth=[std_rhos(2)/max_rhos, std_muKs(2)/max_muKs, 0];  %Varability in c, V, K
                sigma_cost=0;  %Variability in cost
                meanCost=.95;

                this_params=randParams(subStrains, T, d, S0, meanCost, cstar, Vstar, Kstar, seg_rate, conj_rate, epsilon, sigma_cost, sigma_growth);
            end
            
             %Initial conditions
            Bf0s=(x0/subStrains).*ones(1,subStrains);  %Inital plasmid-bearing type frequency (uniform)
            Bp0s=(x0/subStrains).*ones(1,subStrains); %zeros(1,subStrains); %1e-5.*Bp0s; %zeros(1,numStrains);
            Bp0s(1)=x0/subStrains;
            ic=[S0, Bp0s, Bf0s];  %initial conditions
            
            %Solve ODE
            [times, ys, dyn, t_end, pf]=simulateExtinctionMany(this_params, ic); % solve ode

            pFractions(n)=pf;
            relFitness(n)=log(ys(end,2:subStrains+1))/log(ys(end,subStrains+2:end));
        end

        disp(['  mean(pFractions)=',num2str(mean(pFractions))]);
        disp(['  std(pFractions)=',num2str(std(pFractions))]);

        expe_pFraction(icr, :)=[mean(pFractions) std(pFractions)];
        expe_relFitness(icr, :)=[mean(relFitness) std(relFitness)];
    end

    %Save file
    if(si==1)    % SAMPLE RANDOM STRAINS
        save([runDir,'conjRate_assay_sigmaData.mat'],'expe_relFitness','expe_pFraction','conj_rates','numExperiments');
    else
        save([runDir,'conjRate_assay_sigma0.mat'],'expe_relFitness','expe_pFraction','conj_rates','numExperiments'); 
    end
    
end
