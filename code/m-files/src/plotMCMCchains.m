
%% PLOT MCMC CHAINS

this_nclones=1;
priors={'uniform'}; %,'lognormal','beta','gamma'

this_min_muKs=max_muKs;
this_min_rhos=max_rhos;
this_max_muKs=0;
this_max_rhos=0;

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
        
        for p=1:length(priors)
            lbl_prior=priors{p};
            
            %Open output files from MCMC and load posterior distributions
            dataDir=[rootDir,lbl_strain,'_',lbl_prior,'_',lbl_clones,'_',lbl_type,'/'];
            
            if exist([dataDir, 'posterior_samples_rho.csv'],'file')
                this_posterior_samples_rho = table2array(readtable([dataDir, 'posterior_samples_rho.csv']));
                this_posterior_samples_muK = table2array(readtable([dataDir, 'posterior_samples_muK.csv']));
                
                this_mean_rho=mean(this_posterior_samples_rho);
                this_mean_muK=mean(this_posterior_samples_muK);
                
                this_std_rho=std(this_posterior_samples_rho);
                this_std_muK=std(this_posterior_samples_muK);
                
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
                
                
                %Plot chains rho
                subaxis(4, 4, 1, 1, 2, 1, 'SpacingVert',0.02,'PaddingRight',0.05); hold all;
                plot(1:length(this_posterior_samples_rho), this_posterior_samples_rho,'-','Color',this_color); hold on;
                xlim([1, length(this_posterior_samples_rho)]);
                ylim([this_min_rhos*.9 this_max_rhos*1.1]);
                xticklabels([]);
                set(gca,'fontsize',20);
                ylabel('\rho','FontSize',24);

                
                %Plot chains muK
                subaxis(4, 4, 1, 2, 2, 1,'PaddingBottom',0.01,'PaddingRight',0.05); hold on;
                plot(1:length(this_posterior_samples_muK), this_posterior_samples_muK,'-','Color',this_color); hold on;
                xlim([1, length(this_posterior_samples_muK)]);
                ylim([this_min_muKs*.9 this_max_muKs*1.1]);
                set(gca,'fontsize',20);
                ylabel('V_{max}/K_m','FontSize',24);
                xlabel('Iteration');
                
            else
                disp(['Not found ',[dataDir, 'posterior_samples_rho.csv']]);
            end
        end
    end
end
