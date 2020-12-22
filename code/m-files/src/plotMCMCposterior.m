
%% PLOT MCMC POSTERIOR DISTRIBUTIONS

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
                
                %Plot 2D distribution
                this_scatter=scatter(this_posterior_samples_muK, this_posterior_samples_rho,'MarkerFaceColor',this_color,'MarkerEdgeColor',this_color); hold on;
                this_scatter.MarkerFaceAlpha = .2;
                this_scatter.MarkerEdgeAlpha = .4;
                
                %TEMP: ANNOTATE SAMPLES
                iMCMC=find(strcmp(MCMC_strains,lbl_strain ));
                sel_muKs=MCMC_muKs(iMCMC);
                sel_rhos=MCMC_rhos(iMCMC);
                plot(sel_muKs(1), sel_rhos(1), 'ok', 'MarkerFaceColor', color_WT); hold on;
                plot(sel_muKs(2), sel_rhos(2), 'ok', 'MarkerFaceColor', color_TC); hold on;
                
                xlim([this_min_muKs*.8 this_max_muKs*1.2]);
                ylim([this_min_rhos*.99 this_max_rhos*1.01]);
                set(gca,'fontsize',20);
                ylabel('\rho','FontSize',24);
                xlabel('V_{max}/K_m','FontSize',24);
                
            else
                disp(['Not found ',[dataDir, 'posterior_samples_rho.csv']]);
            end
        end
    end
    
end

