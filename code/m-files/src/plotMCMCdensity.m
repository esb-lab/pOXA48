

%% SIMULATE AND PLOT PAIR-WISE COMPARISON BETWEEN TC AND WT (MCMC PARAMETERS)

num_samples=100;
nclones=1;

%Load plasmid-free data
fileName_WT=['ODdata_',lbl_strain,'_WT.csv'];
this_data_WT = csvread([strainPath, fileName_WT],1);
this_times_WT=this_data_WT(:,1);
this_meanOD_WT=this_data_WT(:,2);
this_stdOD_WT=this_data_WT(:,3);
this_meanCFUs_WT=OD2CFU(this_meanOD_WT);
this_stdCFUs_WT=OD2CFU(this_stdOD_WT);

%Load plasmid-bearer data
fileName_TC=['ODdata_',lbl_strain,'_TC.csv'];
this_data_TC = csvread([strainPath, fileName_TC],1);
this_times_TC=this_data_TC(:,1);
this_meanOD_TC=this_data_TC(:,2);
this_stdOD_TC=this_data_TC(:,3);
this_meanCFUs_TC=OD2CFU(this_meanOD_TC);
this_stdCFUs_TC=OD2CFU(this_stdOD_TC);

% Random sample of MCMC posterior distribution
for t=1:length(all_plasmids)
    istrain=find(strcmp(all_strains,selected_strains{s}));
    if ~isempty(istrain)
        code_strain=all_codes{istrain};
        lbl_strain=selected_strains{s};
        lbl_type=all_plasmids{t};
        lbl_prior='uniform';
        lbl_clones=['nclones',num2str(this_nclones)];
        
        
        if t==1
            this_color=color_WT;
            x0=this_meanCFUs_WT(1);
        else
            this_color=color_TC;
            x0=this_meanCFUs_TC(1);
        end
        
        this_params.strains{t}='mcmc';
        this_params.species{t}='sim';
        if t==1
            this_params.plasmids{s}='TC';
            
        else
            this_params.plasmids{s}='WT';
        end
        
        %Open output files from MCMC and load posterior distributions
        dataDir=[rootDir,lbl_strain,'_',lbl_prior,'_',lbl_clones,'_',lbl_type,'/'];
        
        if exist([dataDir, 'posterior_samples_rho.csv'],'file')
            %disp(['Loading ',dataDir, 'posterior_samples_rho.csv']);
            this_posterior_samples_rho = table2array(readtable([dataDir, 'posterior_samples_rho.csv']));
            this_posterior_samples_muK = table2array(readtable([dataDir, 'posterior_samples_muK.csv']));
            
            
            for ns=1:num_samples
                %Sample posterior distribution
                isample=randi(length(this_posterior_samples_rho),1);
                this_sample_rho=this_posterior_samples_rho(isample);
                this_sample_muK=this_posterior_samples_muK(isample);
                
                this_params.seg_rate=seg_rate;
                this_params.conj_rate=conj_rate;
                this_params.S0=S0;
                this_params.T=T;
                this_params.d=d;
                this_params.epsilon=epsilon;
                this_params.numStrains=1;
                this_params.cs=[this_sample_rho 0];
                this_params.Vs=[this_sample_muK 0];
                this_params.Ks=[1 1];
                
                %disp([num2str(isample),':     rho=',num2str(this_sample_rho),' muK=',num2str(this_sample_muK)])
                
                ic=[S0, x0, 0];  %initial conditions
                
                %Solve ODE
                options.RelTol = 1e-12;
                options.AbsTol = 1e-12;
                [this_times, this_ys] = ode15s(@(t,x)fMany(t,x, this_params),[0,1],ic, options);
                
                %Plot realizations
                S=this_ys(:,1);
                B=this_ys(:,2:end);
                maxY=1e9;
                %maxY=max([max(sum(B(:,1:2),2)), max(sum(B(:,2:end),2)) ]);
                plot(this_times*T, B(:,1),'-', 'Color', this_color,'LineWidth',1); hold on;
                ylim([0 14e8]);
                alpha(0.2);
                ylabel('Density','fontsize',24);
                xlim([0 this_times_TC(end)]);
                xlabel('Time','fontsize',24);
                if T==24
                    xticks(0:T/6:T);
                end
                set(gca,'FontSize',20);
            end
        end
    end
end


