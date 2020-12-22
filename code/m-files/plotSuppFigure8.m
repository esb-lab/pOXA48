% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Supplementary Figure 8 in
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

seg_rate=1e-8;  %Segregation rate
conj_rate=0e-12;  %Conjugation rate

epsilon=1e-6; %extinction threshold

max_muKs=10e-10; %Limits for gmmodel
max_rhos=12e8;

color_TC=[0     0.56863     0.56863];
color_WT=[0.83137     0.36863           0];


%% SIMULATION PARAMETERS

strains={ 'C001',  'C002',  'C006',  'C011',  'C012',  'C021',  'C022',  'C031',  'C051',  'C063',  'C094',  'C107',  'C115',  'C131',  'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',  'K141',  'K168',  'K177',  'K200',  'K201',  'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
plasmids={'WT','TC'};

totStrains=length(strains);  % Number of strains (data)
simStrains=1e3;  %Number of strains (sample gmmodel)

%For conjugation assay
numExperiments=1e4;
subStrains=2; %Size of sub-community

%% PATHS AND FILE STRUCTURE


mcmcPath='../../data/MCMC_params.mat';
figDir='../../figures/';


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

%Colormap
cmap_brewer=cbrewer('qual', 'Paired', totStrains*2);
cmap=cmap_brewer([1:2:totStrains*2 2:2:totStrains*2], :);

%Load parameters
mcmc_params=dataParameters(mcmcPath, strains, cmap, T, d, S0, seg_rate, conj_rate, epsilon);


%% SIMULATE AND PLOT PAIR-WISE COMPARISON BETWEEN TC AND WT (DATA PARAMETERS)

disp([newline,'====== Supp Figure 8 ======',newline]);

relFitness=zeros(1,totStrains);
for istrain=1:totStrains
    
    fprintf('.');
    if mod(istrain,totStrains)==0
        disp([num2str(istrain),' ']);
    end
    
    Bp0s=x0/2;  %Inital plasmid-bearing type frequency (uniform)
    Bf0s=x0/2; %Initial plasmid-bearing population
    
    ic=[S0, Bp0s, Bf0s];  %initial conditions
    
    %Parameter values
    this_strains={strains{istrain}};
    strain_params=dataParameters(mcmcPath, this_strains, [cmap(istrain,:)], T, d, S0, seg_rate, conj_rate, epsilon);
    
    if ~strcmp(strain_params.strains{1},'X')
        [times, ys, dyn, t_end, pf]=simulateExtinctionMany(strain_params, ic); % solve ode
        relFitness(istrain)=log(ys(end,2))/log(ys(end,3)); %pair-wise comparison
        
    else
        strain_params;
    end
    
end

%% Compare RelFitness experiment vs model



figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);

plot([0 2], [1 1], ':k','LineWidth',1); hold on;
Ar=[];
Br=[];
for istrain=1:totStrains
    idata=find(strcmp(DATA_strains, strains{istrain}));
    if idata>0
        this_color=cmap(istrain,:);
        plot(DATA_relFitness(idata), relFitness(istrain), 'o', 'MarkerFaceColor', this_color, 'MarkerEdgeColor', this_color); hold on;
        
        Ar=[Ar DATA_relFitness(idata)];
        Br=[Br relFitness(istrain)];
    end
    
end

b = polyfit(Ar, Br, 1);
f = polyval(b, Ar);
Bbar = mean(Br);
SStot = sum((Br - Bbar).^2);
SSreg = sum((f - Bbar).^2);
SSres = sum((Br - f).^2);
R2 = 1 - SSres/SStot;
R = corrcoef(Ar,Br);
Rsq = R(1,2).^2;

y1=.75;
y2=1.25;
axis([y1 y2 y1 y2]);

xq = linspace(y1,y2,100);
yq = polyval(b,xq);
plot(xq,yq, 'k--','LineWidth',1);

text(xq(end),yq(end),['R^2=',num2str(Rsq)],'FontSize',20,'HorizontalAlignment','Right','VerticalAlignment','Bottom');

xticks(linspace(y1, y2, 3));
yticks(linspace(y1, y2, 3));

set(gca,'FontSize',20);
xlabel('Relative fitness (data)', 'FontSize', 24);
ylabel('Relative fitness (model)', 'FontSize', 24);

eval(['export_fig ',figDir,'SuppFigure8.pdf']);


