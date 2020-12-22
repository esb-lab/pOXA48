% This is the repository for the Matlab codes of the numerical simulations
% of plasmid dynamics in complex communities. The scripts can be used to
% generate Figure 5a and Figure 5b of the article
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

%Fit bivariate Normal distribution to TC/WT parameters
gm_TC=getGMModel(MCMC_muKs(iTC), MCMC_rhos(iTC), max_muKs, max_rhos);
gm_WT=getGMModel(MCMC_muKs(iWT), MCMC_rhos(iWT), max_muKs, max_rhos);

mean_muKs=[mean(MCMC_muKs(iTC)) mean(MCMC_muKs(iWT))];
mean_rhos=[mean(MCMC_rhos(iTC)) mean(MCMC_rhos(iWT))];
disp(['Mean muKs (TC,WT)=(',num2str(mean_muKs(1)),',',num2str(mean_muKs(2)),')']);
disp(['Mean rhos (TC,WT)=(',num2str(mean_rhos(1)),',',num2str(mean_rhos(2)),')']);

%Colormap
cmap_brewer=cbrewer('qual', 'Paired', totStrains*2);
cmap=cmap_brewer([1:2:totStrains*2 2:2:totStrains*2], :);

%Load parameters
mcmc_params=dataParameters(mcmcPath, strains, cmap, T, d, S0, seg_rate, conj_rate, epsilon);


%% SAMPLE RANDOM STRAINS

gm_sample_WT=(random(gm_WT,simStrains)).*[max_muKs, max_rhos];
gm_sample_TC=(random(gm_TC,simStrains)).*[max_muKs, max_rhos];
gm_sample=[gm_sample_TC; gm_sample_WT];

gm_params=gmParameters(gm_sample,  T, d, S0, seg_rate, conj_rate, epsilon);

[h1_muK,p1_muK]=ttest2(gm_sample_TC(:,1), gm_sample_WT(:,1), .05, 'both','unequal');
if h1_muK==1
    disp([newline,'There are differences in muK between TC and WT distributions (p-value=',num2str(p1_muK),')']);
else
    disp([newline,'There are no significant differencesin muK between TC and WT distributions (p-value=',num2str(p1_muK),')']);
end

[h1_rho,p1_rho]=ttest2(gm_sample_TC(:,2), gm_sample_WT(:,2), .05, 'both','unequal');
if h1_rho==1
    disp([newline,'There are differences in rho between TC and WT distributions (p-value=',num2str(p1_rho),')', newline]);
else
    disp([newline,'There are no significant differences in rho between TC and WT distributions (p-value=',num2str(p1_rho),')', newline]);
end
%% PLOT 2-D PARAMETER DISTRIBUTION


disp([newline,'====== Figure 5a ======',newline]);

figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);

plotParams(mcmc_params, max_muKs, max_rhos);

%Contour plots
gm_TC=getGMModel(MCMC_muKs(iTC), MCMC_rhos(iTC), max_muKs, max_rhos);
gm_WT=getGMModel(MCMC_muKs(iWT), MCMC_rhos(iWT), max_muKs, max_rhos);

fc_TC=fcontour(@(x,y)reshape(pdf(gm_TC,[x(:),y(:)]),size(x)),[0 1],'LineWidth',1,'LineColor',color_TC); hold on;
fc_WT=fcontour(@(x,y)reshape(pdf(gm_WT,[x(:),y(:)]),size(x)),[0 1],'LineWidth',1,'LineColor',color_WT,'LineStyle','-'); hold on;

fc_TC.LineWidth = 2;
fc_TC.LineStyle = '-';
fc_TC.LevelList = [1];

fc_WT.LineWidth = 2;
fc_WT.LineStyle = '-';
fc_WT.LevelList = [1];

nticks=5;

xticks(0:1/nticks:1);
yticks(0:1/nticks:1);

mean_muKs=[mean(MCMC_muKs(iTC)) mean(MCMC_muKs(iWT))];
mean_rhos=[mean(MCMC_rhos(iTC)) mean(MCMC_rhos(iWT))];
disp(['Mean muKs (TC,WT)=(',num2str(mean_muKs(1)),',',num2str(mean_muKs(2)),')']);
disp(['Mean rhos (TC,WT)=(',num2str(mean_rhos(1)),',',num2str(mean_rhos(2)),')']);
axis([.1 1. .4 1]);

eval(['export_fig ',figDir,'Figure5a.pdf']);

%% SIMULATE AND PLOT PAIR-WISE COMPARISON BETWEEN TC AND WT (DATA PARAMETERS)

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


%% SIMULATE AND PLOT PAIR-WISE COMPARISON BETWEEN TC AND WT (SAMPLE PARAMS)

gm_relFitness=zeros(1,simStrains);
for istrain=1:simStrains
    
    fprintf('.');
    if mod(istrain,1e2)==0
        disp([num2str(istrain),' ']);
    end
    
    Bp0s=x0/2;  %Inital plasmid-bearing type frequency (uniform)
    Bf0s=x0/2; %Initial plasmid-bearing population
    
    ic=[S0, Bp0s, Bf0s];  %initial conditions
    
    %Parameter values
    sub_strains={gm_params.strains{istrain}};
    pair_params=subParameters(gm_params, sub_strains);
    
    if ~strcmp(pair_params.strains{1},'X')
        [times, ys, dyn, t_end, pf]=simulateExtinctionMany(pair_params, ic); % solve ode
        gm_relFitness(istrain)=log(ys(end,2))/log(ys(end,3)); %pair-wise comparison
    end
    
end

%% Plot DFE (DATA PARAMS & SAMPLE PARAMS)


disp([newline,'====== Figure 5b ======',newline]);

figure(); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'FontSize',24);

bin_min=0.5;
bin_max=1.5;
bin_delta_sample=bin_max/51;
ebins_sample=bin_min:bin_delta_sample:bin_max;

%SAMPLE DIST
h_sample=histogram(gm_relFitness,ebins_sample,  'FaceColor','y', 'FaceAlpha',1,'EdgeColor','w');
sample_y=h_sample.Values;
sample_x=ebins_sample(2:end)-bin_delta_sample/2;
disp(['Mean(gm_relFitness)=',num2str(mean(gm_relFitness))]);
disp(['Std(gm_relFitness)=',num2str(std(gm_relFitness))]);

%DATA DIST
bin_delta_data=bin_max/51;
ebins_data=bin_min:bin_delta_data:bin_max;
h_data=histogram(relFitness,ebins_data,  'FaceColor','k', 'FaceAlpha',1,'EdgeColor','w');
data_y=h_data.Values;
data_x=ebins_data(2:end)-bin_delta_data/2;
mean_relFitness=mean(relFitness);
std_relFitness=std(relFitness);
disp(['mean(relFitness)=',num2str(mean_relFitness)]);
disp(['std(relFitness)=',num2str(std_relFitness)]);

bar(data_x, data_y/max(data_y) ,'LineWidth',1,'FaceColor',[.9 .9 .9]); hold on;
plot(sample_x, sample_y/max(sample_y) ,'k-','LineWidth',2);

plot([1 1],[0 sample_x(end)*1.1],'-k','LineWidth',1);
xlim([0.5 1.5]);
ylim([0 1]);
xticks(0:.2:2);

set(gca,'FontSize',20);
xlabel('Relative fitness','FontSize',24);
ylabel('Frequency','FontSize',24);

eval(['export_fig ',figDir,'Figure5b.pdf']);


