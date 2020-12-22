clc
close all
clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');

%%
figDir='../../figures/src/WTvTC/';
strainPath='../../data/ODdata_strains/';
dataPath='../../data/MCMC_params.mat';

%% CREATE DIRECTORY STRUCTURE

if ~isempty(figDir)
    if ~exist(figDir,'dir')
        mkdir(figDir);
        disp(['mkdir ',figDir]);
    end
end

%% PARAMETERS

T=24;   %Length of experiment
S0=1;    %Resource concentration
x0=1e6;  %Initial bacterial density
d=0.0;  %dilution rate
seg_rate=1e-8;  %Segregation rate
conj_rate=0e-12;  %Conjugation rate
epsilon=1e-6; %extinction threshold

codes={ 'Ec01',  'Ec02',  'Ec03',  'Ec04',  'Ec05',  'Ec06',  'Ec07',  'Ec08',  'Ec09',  'Ec10',  'Ec11',  'Ec12',  'Ec13',  'Ec14',  'Ec15',  'Ec16',  'Ec17',  'Ec18',  'Ec19',  'Ec20',  'Ec21',  'Ec22',  'Ec23',  'Ec24',  'Ec25',  'Kpn01',  'Kpn02',  'Kpn03',  'Kpn04',  'Kpn05',  'Kpn06',  'Kpn07',  'Kpn08',  'Kpn09',  'Kpn10',  'Kpn11',  'Kpn12',  'Kpn13',  'Kpn14',  'Kpn15',  'Kpn16',  'Kpn17',  'Kpn18',  'Kpn19',  'Kpn20',  'Kpn21',  'Kpn22',  'Kpn23',  'Kpn24',  'Kpn25'};
strains={ 'C001',  'C002',  'C006',  'C011',  'C012',  'C021',  'C022',  'C031',  'C051',  'C063',  'C094',  'C107',  'C115',  'C131',  'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',  'K141',  'K168',  'K177',  'K200',  'K201',  'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};


%strains={'K177', 'C001', 'C022', 'C002',  'C006',  'C011',  'C012',  'C021',  'C031',  'C051',  'C063',    'C107',  'C115',   'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',    'K168',  'K200',    'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
plasmids={'WT','TC'};
load(dataPath);
disp([num2str(length(MCMC_strains)),' MCMC output files loaded']);

max_muKs=10e-10;
max_rhos=12e8;

totStrains=length(strains);

cmap_brewer=cbrewer('qual', 'Paired', totStrains*2);
cmap=cmap_brewer([1:2:totStrains*2 2:2:totStrains*2], :);

iTC=find(strcmp(MCMC_plasmids,'TC' ));
iWT=find(strcmp(MCMC_plasmids,'WT' ));

iK=find(strcmp(MCMC_species,'K' ));
iE=find(strcmp(MCMC_species,'E' ));

data_params=dataParameters(dataPath, strains, cmap, T, d, S0, seg_rate, conj_rate, epsilon);


%%

for istrain=1:totStrains
    this_strain=strains{istrain};
    imcmc=find(strcmp(MCMC_strains,strains{istrain} ));

    muK_TC=MCMC_muKs(imcmc(1));  %muK <- TC
    rho_TC=MCMC_rhos(imcmc(1));  %rho <- TC

    muK_WT=MCMC_muKs(imcmc(2));  %muK <- WT
    rho_WT=MCMC_rhos(imcmc(2));  %rho <- WT

    data_params.Vs(istrain)=muK_TC;
    data_params.cs(istrain)=rho_TC;
    data_params.Vs(istrain+totStrains)=muK_WT;
    data_params.cs(istrain+totStrains)=rho_WT;

    disp('_______');
    disp([num2str(istrain),': ',data_params.strains{istrain}]);
    disp(['  muK (TC,WT)=(',num2str(data_params.Vs(istrain)),',',num2str(data_params.Vs(istrain+totStrains)),')']);
    disp(['  rho (TC,WT)=(',num2str(data_params.cs(istrain)),',',num2str(data_params.cs(istrain+totStrains)),')']);

    %Load plasmid-free data
    fileName_WT=['ODdata_',this_strain,'_WT.csv'];
    this_data_WT = csvread([strainPath, fileName_WT],1);
    this_times_WT=this_data_WT(:,1);
    this_meanOD_WT=this_data_WT(:,2);
    this_stdOD_WT=this_data_WT(:,3);
    this_meanCFUs_WT=OD2CFU(this_meanOD_WT);
    this_stdCFUs_WT=OD2CFU(this_stdOD_WT);

    %Load plasmid-bearer data
    fileName_TC=['ODdata_',this_strain,'_TC.csv'];
    this_data_TC = csvread([strainPath, fileName_TC],1);
    this_times_TC=this_data_TC(:,1);
    this_meanOD_TC=this_data_TC(:,2);
    this_stdOD_TC=this_data_TC(:,3);
    this_meanCFUs_TC=OD2CFU(this_meanOD_TC);
    this_stdCFUs_TC=OD2CFU(this_stdOD_TC);

    %Plot data
    figure(1); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'fontsize',18);
    set(gcf,'Units','normalized','Position',[0. 0.25 .2 .3])
    p1_WT=plot(this_times_WT, this_meanCFUs_WT,'k-'); hold on;  %Plasmid-free
    plot(this_times_WT, this_meanCFUs_WT+this_stdCFUs_WT,'k:','LineWidth',1); hold on;  %Plasmid-free
    plot(this_times_WT, this_meanCFUs_WT-this_stdCFUs_WT,'k:','LineWidth',1); hold on;  %Plasmid-free
    
    p1_TC=plot(this_times_TC, this_meanCFUs_TC,'r-'); hold on;
    plot(this_times_TC, this_meanCFUs_TC+this_stdCFUs_TC,'r:','LineWidth',1); hold on;  %Plasmid-free
    plot(this_times_TC, this_meanCFUs_TC-this_stdCFUs_TC,'r:','LineWidth',1); hold on;  %Plasmid-free
    
    set(gca,'fontsize',18);
    ylabel('Density (cells/ml)','fontsize',24);
    xlabel('Time (hours)','fontsize',24);
    xticks(0:6:24);
    xlim([0, 24])
    title([this_strain,' (data)'],'fontsize',18);
    ylim([0 12e8]);
    legend([p1_WT, p1_TC],{'Plasmid-free','Plasmid-bearing'},'Location','SouthEast');

    
    %Now import parameters from MCMC_params.mat
    iTC=find(strcmp(MCMC_plasmids,'TC' ));
    iWT=find(strcmp(MCMC_plasmids,'WT' ));
    imcmc=find(strcmp(MCMC_strains,this_strain ));

    istrain_TC=intersect(iTC, imcmc);
    istrain_WT=intersect(iWT, imcmc);

    disp(['_____ ',this_strain]);
    disp(['WT:']);
    disp(['** mean(muK)=',num2str(muK_WT)]);
    disp(['** mean(rho)=',num2str(rho_WT)]);

    disp(['TC:']);
    disp(['** mean(muK)=',num2str(muK_TC)]);
    disp(['** mean(rho)=',num2str(rho_TC)]);

    params.seg_rate=seg_rate;
    params.conj_rate=conj_rate;
    params.S0=S0;
    params.T=T;
    params.d=d;
    params.epsilon=epsilon;
    params.Vs=[muK_TC muK_WT];
    params.Ks=[1 1];
    params.cs=[rho_TC rho_WT];
    params.numStrains=1;

    %Plot growth curves (model)
    figure(2); clf('reset'); set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white'); set(gca,'fontsize',18);
    set(gcf,'Units','normalized','Position',[0. 0.25 .2 .3])

    ic_WT=[1, 0, this_meanCFUs_WT(1)];
    options = odeset('NonNegative', 1:3);
    [sim_times_WT, sim_y_WT] = ode23(@(t,x)fMany(t,x, params),[0,1],ic_WT,options);
    p_WT=plot(sim_times_WT.*params.T, sim_y_WT(:,3),'r-'); hold on; %Plasmid-free
    
    ic_TC=[1, this_meanCFUs_TC(1), 0];
    [sim_times_TC, sim_y_TC] = ode23(@(t,x)fMany(t,x, params),[0,1],ic_TC,options);
    p_TC=plot(sim_times_TC.*params.T, sim_y_TC(:,2),'k-'); hold on; %Plasmid-bearing
    set(gca,'fontsize',18);
    legend([p_WT, p_TC],{'Plasmid-free','Plasmid-bearing'},'Location','SouthEast');
    title([this_strain,' (model)'],'fontsize',18);
    ylabel('Density (cells/ml)','fontsize',24);
    xlabel('Time (hours)','fontsize',24);
    xticks(0:6:24);
    xlim([0, 24])
    ylim([0 12e8]);

    figure(1);
    eval(['export_fig ',figDir,'',this_strain,'_OD600_data.pdf']);
    close;
    
    figure(2);
    eval(['export_fig ',figDir,'',this_strain,'_OD600_model.pdf']);
    close;

end

