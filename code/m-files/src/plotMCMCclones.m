
%% DATA CLONING


priors={'uniform'};
nclones=[1, 5, 10, 20];

for t=1:length(all_plasmids)
    
    if t==1
        this_color=color_WT;
    else
        this_color=color_TC;
    end
    istrain=find(strcmp(all_strains,selected_strains{s}));
    if ~isempty(istrain)
        code_strain=all_codes{istrain};
        lbl_strain=selected_strains{s};
        lbl_type=all_plasmids{t};
        
        
        mean_nclones=zeros(length(nclones), 2);
        std_nclones=zeros(length(nclones), 2);
        
        N=5e4;
        M_rho=zeros(N, length(nclones));
        M_muK=zeros(N, length(nclones));
        
        for n=1:length(nclones)
            
            lbl_clones=['nclones',num2str(nclones(n)),'_'];
            
            for p=1:length(priors)
                lbl_prior=priors{p};
                
                
                %Open output files from MCMC and load posterior distributions
                dataDir=[rootDir,lbl_strain,'_',lbl_prior,'_',lbl_clones,'',lbl_type,'/'];
                
                if exist([dataDir, 'posterior_samples_rho.csv'],'file')
                    
                    this_posterior_samples_rho = table2array(readtable([dataDir, 'posterior_samples_rho.csv']));
                    this_posterior_samples_muK = table2array(readtable([dataDir, 'posterior_samples_muK.csv']));
                    
                    
                    this_mean_rho=mean(this_posterior_samples_rho);
                    this_mean_muK=mean(this_posterior_samples_muK);
                    this_std_rho=std(this_posterior_samples_rho);
                    this_std_muK=std(this_posterior_samples_muK);
                    
                    mean_nclones(n,:)=[this_mean_rho this_mean_muK];
                    std_nclones(n,:)=[this_std_rho this_std_muK];
                    
                    M_muK(1:N, n)=this_posterior_samples_muK(end-N+1:end);
                    M_rho(1:N, n)=this_posterior_samples_rho(end-N+1:end);
                    
                else
                    disp(['Not found ',dataDir, '']);
                end
            end
        end
        
        
        if sum(sum(mean_nclones))>0
            
            %Plot rho
            subaxis(4, 4, 4, 4,1,1,'PaddingTop',0.0,'PaddingBottom',0.0,'SpacingHoriz',0.05);
            
            plot(nclones(std_nclones(:,1)>0),std_nclones(std_nclones(:,1)>0,1),'o-','Color',this_color,'MarkerFaceColor',this_color); hold on;
            
            %ylim([mean(mean_nclones(mean_nclones(:,1)>0,1))-max(std_nclones(:,1))*3 mean(mean_nclones(mean_nclones(:,1)>0,1))+max(std_nclones(:,1))*3]);
            
            ylim([0 this_max_rhos]);
            %xticks(1:length(nclones));
            %xlim([0. length(nclones)+.5]);
            xlim([0. 1.1*max(nclones)]);
            xticks(nclones);
            set(gca,'fontsize',20);
            axis tight
            %xticklabels(nclones)
            ylabel('std(\rho)','FontSize',24);
            xlabel('Number of clones');
            
            %Plot muK
            subaxis(4, 4, 3, 4,1,1,'PaddingTop',0.0,'PaddingBottom',0.0,'SpacingHoriz',0.05);
            plot(nclones(std_nclones(:,2)>0),std_nclones(std_nclones(:,2)>0,2),'o-','Color',this_color,'MarkerFaceColor',this_color); hold on;
            
            axis tight
            xticks(nclones);
            xlim([0. 1.1*max(nclones)]);
            set(gca,'fontsize',20);
            %xticklabels(nclones)
            ylabel('std(V_{max}/K_m)','FontSize',24);
            xlabel('Number of clones');
            
        end
    end
end

%% ANOVA

[p_muK,~] = anova1(M_muK,[],'off');
[p_rho,~] = anova1(M_rho,[],'off');

if p_muK>.005 && p_rho>0.005
    disp([code_strain, ': No significant differences between different number of clones (ANOVA)']);
else
    disp([code_strain, ': Significant differences between different numbers of clones (ANOVA)']);
end

disp(['   muK: p-value=',num2str(p_muK)]);
disp(['   rho: p-value=',num2str(p_rho)]);


