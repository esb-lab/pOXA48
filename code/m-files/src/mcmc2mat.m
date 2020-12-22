clc
close all
clear all

run('lib/addpath_recurse')
addpath_recurse('src/');
addpath_recurse('lib/');



%%


codes={ 'Ec01',  'Ec02',  'Ec03',  'Ec04',  'Ec05',  'Ec06',  'Ec07',  'Ec08',  'Ec09',  'Ec10',  'Ec11',  'Ec12',  'Ec13',  'Ec14',  'Ec15',  'Ec16',  'Ec17',  'Ec18',  'Ec19',  'Ec20',  'Ec21',  'Ec22',  'Ec23',  'Ec24',  'Ec25',  'Kpn01',  'Kpn02',  'Kpn03',  'Kpn04',  'Kpn05',  'Kpn06',  'Kpn07',  'Kpn08',  'Kpn09',  'Kpn10',  'Kpn11',  'Kpn12',  'Kpn13',  'Kpn14',  'Kpn15',  'Kpn16',  'Kpn17',  'Kpn18',  'Kpn19',  'Kpn20',  'Kpn21',  'Kpn22',  'Kpn23',  'Kpn24',  'Kpn25'};
strains={ 'C001',  'C002',  'C006',  'C011',  'C012',  'C021',  'C022',  'C031',  'C051',  'C063',  'C094',  'C107',  'C115',  'C131',  'C141',  'C201',  'C227',  'C232',  'C247',  'C261',  'C286',  'C290',  'C302',  'C309',  'C324',  'K037',  'K038',  'K087',  'K094',  'K112',  'K114',  'K125',  'K141',  'K168',  'K177',  'K200',  'K201',  'K209',  'K213',  'K216',  'K224',  'K225',  'K241',  'K248',  'K249',  'K253',  'K257',  'K275',  'K285',  'K300'};
plasmids={'WT','TC'};

ucode=uniqueStrCell(codes);

dataDir='../../data/runs_mcmc/strains/';
matDir='../../data/';
figDir='../../figures/src/';

%%
clc
str_muK='mean(muK)=';
str_rho='mean(rho)=';

MCMC_species={};
MCMC_strains={};
MCMC_plasmids={};
MCMC_muKs=[];
MCMC_rhos=[];
for istrain=1:length(strains)
    
    
    for iplasmid=1:length(plasmids)

        logPath_txt=[dataDir,'',strains{istrain},'_',plasmids{iplasmid},'/log.txt'];
        if exist(logPath_txt,'file')
            logPath=logPath_txt;
        end
        
        if exist(logPath,'file')
            disp([num2str(istrain),': Loading ',logPath]);

            fid=fopen(logPath);
            tline = fgetl(fid);
            val_muK=[];
            val_rho=[];
            while ischar(tline)
                if contains(tline, str_muK)==1
                    val_muK=str2double(extractAfter(tline,str_muK));
                    MCMC_muKs=[MCMC_muKs, val_muK];
                end

                if contains(tline, str_rho)==1
                    val_rho=str2double(extractAfter(tline,str_rho));
                    MCMC_rhos=[MCMC_rhos, val_rho];

                end

                tline = fgetl(fid);
            end
            fclose(fid);


            if ~isempty(val_muK) 
                MCMC_strains{length(MCMC_muKs)}=strains{istrain};
                MCMC_plasmids{length(MCMC_muKs)}=plasmids{iplasmid};
                if strfind(strains{istrain},'K')
                    MCMC_species{length(MCMC_muKs)}='K';
                else
                    MCMC_species{length(MCMC_muKs)}='E';
                end
            else
                disp([num2str(istrain),': *Not found ',logPath]);
                    
            end
        else
            disp([num2str(istrain),': Not found ',logPath]);
        end

    end
    
end

iTC=find(strcmp(MCMC_plasmids,'TC' ));
iWT=find(strcmp(MCMC_plasmids,'WT' ));

iK=find(strcmp(MCMC_species,'K' ));
iE=find(strcmp(MCMC_species,'E' ));

disp([num2str(length(MCMC_strains)),' MCMC parameters loaded']);

%% TODO: Check outliers!



%max_muKs=max(MCMC_muKs);
%max_rhos=max(MCMC_rhos);

max_muKs=10e-10;
max_rhos=12e8;


iok=intersect(find(MCMC_muKs<max_muKs) ,find(MCMC_rhos<max_rhos)) ;

nMCMC_muKs=MCMC_muKs./max_muKs;
nMCMC_rhos=MCMC_rhos./max_rhos;

figure(1); clf('reset');
set(gcf,'DefaultLineLineWidth',2); set(gcf, 'color', 'white');
        
X = [nMCMC_muKs(iok)' nMCMC_rhos(iok)'];
GMModel = fitgmdist(X,1);

nfitMu=GMModel.mu;
nfitSigma=GMModel.Sigma;
fitMu=nfitMu.*[max_muKs max_rhos]
fitSigma=nfitSigma.*[max_muKs max_rhos]

%gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1)); g = gca;
%fcontour(gmPDF,[g.XLim g.YLim])

gm = gmdistribution(nfitMu,nfitSigma);
level_set=[1];
level_color=[.7 .7 .7];
gmPDF=@(x,y)reshape(pdf(gm,[x(:),y(:)]),size(x));
fcontour(gmPDF,[0 1],'LineWidth',.5,'LineColor',level_color); hold on;


color_K=[0 0 1];
color_E=[1 0 0];

iTC_K=intersect(iK, iTC);
p_TC_K=plot(nMCMC_muKs(iTC_K),nMCMC_rhos(iTC_K), 'o','LineWidth',1,'MarkerEdgeColor',color_K); hold on;

iWT_K=intersect(iK, iWT);
p_WT_K=plot(nMCMC_muKs(iWT_K), nMCMC_rhos(iWT_K), 'o','MarkerFaceColor',color_K,'LineWidth',1,'MarkerEdgeColor',color_K); hold on;


iTC_E=intersect(iE, iTC);
p_TC_E=plot(nMCMC_muKs(iTC_E),nMCMC_rhos(iTC_E), 'o','LineWidth',1,'MarkerEdgeColor',color_E); hold on;

iWT_E=intersect(iE, iWT);
p_WT_E=plot(nMCMC_muKs(iWT_E), nMCMC_rhos(iWT_E), 'o','MarkerFaceColor',color_E,'LineWidth',1,'MarkerEdgeColor',color_E); hold on;

nticks=5;
axis([0 1 0 1]);

xticks(0:1/nticks:1);
xticklabels(0:max_muKs/nticks:max_muKs);

yticks(0:1/nticks:1);
yticklabels(0:max_rhos/nticks:max_rhos);


set(gca,'FontSize',16);

legend([p_WT_K, p_TC_K,p_WT_E, p_TC_E],{'Klebsiella','Klebsiella (plasmid-free)', 'Escherichia', 'Escherichia (plasmid-free)'},'FontSize',16, 'Location','SouthWest');
legend boxoff
xlabel('Specific affinity (\mu/K)','FontSize',18);
ylabel('Cell efficiency (\rho)','FontSize',18);

eval(['export_fig ',figDir,'/model_params_mcmc.pdf']);

%%

save([matDir,'MCMC_params.mat'],'MCMC_species','MCMC_strains','MCMC_plasmids','MCMC_muKs','MCMC_rhos')
