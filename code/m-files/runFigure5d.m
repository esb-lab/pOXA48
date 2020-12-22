% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate the data necessary to produce Figure 5d of the manuscript
% "Variability of plasmid fitness effects contributes to plasmid persistence in bacterial communities."
%
% December 16, 2020
% rpm@ccg.unam.mx


clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');


disp([newline,'====== Run Figure 5d ======',newline]);

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


%% SIMULATION PARAMETERS

%From data
strains={'C001', 'C022', 'C002',  'C006',  'C011',  'C012',  'C021',  'C031',  'C051',  'C063',    'C107',  'C115',   'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',    'K168',  'K177', 'K200',    'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
plasmids={'WT','TC'};
totStrains=length(strains);  % Number of strains (data)

%For conjugation assay
numExperiments=1e2; %5e4
subStrains=2; %Size of sub-community

conj_rates=logspace(-11, -8, 101);  %Conjugation rates used to determine optimal 

mus=linspace(0., 0.05, 51); %DFE mean (cost=1-mu)
sigmas=sqrt([0 0.0075 0.015]);  %DFE standard deviation 



%% PATHS AND FILE STRUCTURE

runDir='../../data/runs/figure5d/';
if ~exist(runDir, 'dir')
    mkdir(runDir);
    disp(['mkdir ',runDir]);
end

mcmcPath='../../data/MCMC_params.mat';

save([runDir,'experiment.mat']);

%% LOAD DATA

%MCMC data
load(mcmcPath);
disp([num2str(length(MCMC_strains)),' MCMC output files loaded']);

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

for mi=1:length(mus)

    mean_cost=1-mus(mi);
    
    for si =1:length(sigmas)

        %Variability in cost
        sigma_cost=sigmas(si);   

        disp([newline, num2str(mi),'/', num2str(length(mus)),': mean_cost=',num2str(mean_cost),', sigma_cost=',num2str(sigma_cost), newline]);                    

        expe_pFraction=zeros(length(conj_rates),2);
        expe_relFitness=zeros(length(conj_rates),2);
        
        
        matFile=[runDir,'conjRate_assay_mu',num2str(round(mus(mi)*1000)),'e-3_sigma',num2str(round(sigmas(si)*1e3)),'e-3.mat'];
        
        if ~exist(matFile, 'file')
            for icr=1:length(conj_rates)

                disp(['--> ',num2str(icr),'/',num2str(length(conj_rates)),': conj_rate=',num2str(conj_rates(icr))]);

                conj_rate=conj_rates(icr);

                relFitness=zeros(1,numExperiments);
                pFractions=zeros(1,numExperiments);
                parfor n=1:numExperiments
                       %RANDOM
                        this_colors=[rand rand rand];
                        cstar=mean_rhos(2);
                        Vstar=mean_muKs(2); 
                        Kstar=1;
                        sigma_growth=[std_rhos(2)/max_rhos, std_muKs(2)/max_muKs, 0];  %Varability in c, V, K

                        this_params=randParams(subStrains, T, d, S0, mean_cost, cstar, Vstar, Kstar, seg_rate, conj_rate, epsilon, sigma_cost, sigma_growth);
                  
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

                expe_pFraction(icr, :)=[mean(pFractions) std(pFractions)];
                expe_relFitness(icr, :)=[mean(relFitness) std(relFitness)];
            end
            
            %Save file
            disp(['Saving ',matFile]);
            save(matFile,'expe_relFitness','expe_pFraction','conj_rates','numExperiments'); 
        end
        
    end

end
