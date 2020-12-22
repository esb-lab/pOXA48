%%  PLOT DIFFERENT PRIORS

this_nclones=1;
priors={'uniform','lognormal','beta','gamma'};

this_min_muKs=max_muKs;
this_min_rhos=max_rhos;
this_max_muKs=0;
this_max_rhos=0;

pprior1=[]; pprior2=[];
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
        
        
        lbl_clones=['nclones',num2str(this_nclones)];
        mean_priors=zeros(length(priors), 2);
        std_priors=zeros(length(priors), 2);
        
        N=5e4;
        M_rho=zeros(N, length(priors));
        M_muK=zeros(N, length(priors));
        
        for p=1:length(priors)
            lbl_prior=priors{p};
            
            
            %Open output files from MCMC and load posterior distributions
            dataDir=[rootDir,lbl_strain,'_',lbl_prior,'_',lbl_clones,'_',lbl_type,'/'];
            
            if exist([dataDir, 'posterior_samples_rho.csv'],'file')
                
                %disp(['Loading ',[dataDir, 'posterior_samples_rho.csv']]);
                
                this_posterior_samples_rho = table2array(readtable([dataDir, 'posterior_samples_rho.csv']));
                this_posterior_samples_muK = table2array(readtable([dataDir, 'posterior_samples_muK.csv']));
                
                this_mean_rho=mean(this_posterior_samples_rho);
                this_mean_muK=mean(this_posterior_samples_muK);
                
                this_std_rho=std(this_posterior_samples_rho);
                this_std_muK=std(this_posterior_samples_muK);
                
                mean_priors(p,:)=[this_mean_rho this_mean_muK];
                std_priors(p,:)=[this_std_rho this_std_muK];
                
                
                M_muK(1:N, p)=this_posterior_samples_muK(end-N+1:end);
                M_rho(1:N, p)=this_posterior_samples_rho(end-N+1:end);
                
                %Find boundaries of 2D plot
                if this_min_muKs > min(this_posterior_samples_muK)
                    this_min_muKs= min(this_posterior_samples_muK);
                end
                if this_min_rhos > min(this_posterior_samples_rho)
                    this_min_rhos= min(this_posterior_samples_rho);
                end
                if this_max_muKs < max(this_posterior_samples_muK)
                    this_max_muKs= max(this_posterior_samples_muK);
                end
                if this_max_rhos < max(this_posterior_samples_rho)
                    this_max_rhos= max(this_posterior_samples_rho);
                end
                
                
            else
                
                disp(['Not found ',[dataDir, '']]);
            end
        end
        
        if sum(sum(mean_priors))>0
            
            %Plot rho
            subaxis(4, 4, 4, 1,1,1, 'SpacingVert',0.02,'PaddingTop',0);
            for p=1:length(priors)
                plot(p,mean_priors(p,1),'o','Color',this_color,'MarkerFaceColor',this_color); hold on;
                plot([p p],[mean_priors(p,1)-std_priors(p,1) mean_priors(p,1)+std_priors(p,1)],'-','Color',this_color,'LineWidth',1); hold on;
                plot([p-d p+d],[mean_priors(p,1)+std_priors(p,1) mean_priors(p,1)+std_priors(p,1)],'-','Color',this_color,'LineWidth',1); hold on;
                plot([p-d p+d],[mean_priors(p,1)-std_priors(p,1) mean_priors(p,1)-std_priors(p,1)],'-','Color',this_color,'LineWidth',1); hold on;
            end
            
            %Find boundaries
            if this_min_rhos > min(this_posterior_samples_rho)
                this_min_rhos= min(this_posterior_samples_rho);
            end
            if this_max_rhos < max(this_posterior_samples_rho)
                this_max_rhos= max(this_posterior_samples_rho);
            end
            ylim([this_min_rhos this_max_rhos]);
            
            xlim([0.5 length(priors)+.5]);
            xticks(1:length(priors));
            set(gca,'fontsize',20);
            xtickangle(-30)
            xticklabels(priors)
            ylabel('\rho','FontSize',24);
            
            %Plot muK
            subaxis(4, 4, 3, 1,1,1, 'SpacingVert',0.02);
            for pk=1:length(priors)
                plot(pk,mean_priors(pk,2),'o','Color',this_color,'MarkerFaceColor',this_color); hold on;
                plot([pk pk],[mean_priors(pk,2)-std_priors(pk,2) mean_priors(pk,2)+std_priors(pk,2)],'-','Color',this_color,'LineWidth',1); hold on;
                plot([pk-d pk+d],[mean_priors(pk,2)-std_priors(pk,2) mean_priors(pk,2)-std_priors(pk,2)],'-','Color',this_color,'LineWidth',1); hold on;
                plot([pk-d pk+d],[mean_priors(pk,2)+std_priors(pk,2) mean_priors(pk,2)+std_priors(pk,2)],'-','Color',this_color,'LineWidth',1); hold on;
                
            end
            
            %Find boundaries
            if this_max_muKs < max(this_posterior_samples_muK)
                this_max_muKs= max(this_posterior_samples_muK);
            end
            if this_min_muKs > min(this_posterior_samples_muK)
                this_min_muKs= min(this_posterior_samples_muK);
            end
            ylim([this_min_muKs this_max_muKs]);
            xticks(1:length(priors));
            xlim([0.5 length(priors)+.5]);
            set(gca,'fontsize',20);
            xticklabels(priors)
            xtickangle(-30)
            ylabel('V_{max}/K_m','FontSize',24);
            
            
        end
    end
end

%% ANOVA

[p_muK,~] = anova1(M_muK,[],'off');
[p_rho,~] = anova1(M_rho,[],'off');

if p_muK>.005 && p_rho>0.005
    disp([code_strain, ': No significant differences between different priors (ANOVA)']);
else
    disp([code_strain, ': Significant differences between different priors (ANOVA)']);    
end

disp(['   muK: p-value=',num2str(p_muK)]);
disp(['   rho: p-value=',num2str(p_rho)]);
