
function params=dataParameters(dataPath, this_strains, this_colors, T, d, S0, seg_rate, conj_rate, epsilon)

    numStrains=length(this_strains);

    params.seg_rate=seg_rate;
    params.conj_rate=conj_rate;
    
    load(dataPath);
    iTC=find(strcmp(MCMC_plasmids,'TC' ));
    iWT=find(strcmp(MCMC_plasmids,'WT' ));
    
    params_species={};
    params_cs=zeros(1,2*numStrains);  %rhos
    params_Vs=zeros(1,2*numStrains);  %Vmax
    for icomm=1:length(this_strains)
        
        this_strain=this_strains{icomm};
        istrain=find(strcmp(MCMC_strains,this_strain));
        this_color=this_colors(icomm,:);
        
        %With plasmid
        istrain_TC=intersect(istrain, iTC);
        if ~isempty(istrain_TC)
            
            this_TC_muKs=MCMC_muKs(istrain_TC);
            this_TC_rhos=MCMC_rhos(istrain_TC);
            this_TC_species=MCMC_species{istrain_TC};
            this_TC_strain=[MCMC_strains{istrain_TC},'_',MCMC_plasmids{istrain_TC}];
            
            
            %disp(strcat(this_strain,"_TC: muK=",num2str(this_TC_muKs),", rho=",num2str(this_TC_rhos), ' (',this_TC_species,')'));
            
            params_cs(icomm)=this_TC_rhos;
            params_Vs(icomm)=this_TC_muKs;
            params_species{icomm}=this_TC_species;
            params_strains{icomm}=this_TC_strain;
            params_colors(icomm,:)=this_color;
        else
            params_cs(icomm)=0;
            params_Vs(icomm)=0;
            params_species{icomm}='X';
            params_strains{icomm}='X';
            params_colors(icomm,:)=[.5 .5 .5];
            
            disp(['Not found: ',this_strain,'_TC']);
        end
        
        %Plasmid free
        istrain_WT=intersect(istrain, iWT);
        if ~isempty(istrain_WT)
            
            this_WT_muKs=MCMC_muKs(istrain_WT);
            this_WT_rhos=MCMC_rhos(istrain_WT);
            this_WT_species=MCMC_species{istrain_WT};
            this_WT_strain=[MCMC_strains{istrain_WT},'_',MCMC_plasmids{istrain_WT}];
            
            
            %disp(strcat(this_strain,"_WT: muK=",num2str(this_WT_muKs),", rho=",num2str(this_WT_rhos), ' (',this_WT_species,')'));
            
            params_cs(icomm+numStrains)=this_WT_rhos;
            params_Vs(icomm+numStrains)=this_WT_muKs;
            params_species{icomm+numStrains}=this_WT_species;
            params_strains{icomm+numStrains}=this_WT_strain;
            params_colors(icomm+numStrains,:)=this_color;
        else
            params_cs(icomm+numStrains)=0;
            params_Vs(icomm+numStrains)=0;
            params_species{icomm+numStrains}='X';
            params_strains{icomm+numStrains}='X';
            params_colors(icomm+numStrains,:)=[.5 .5 .5];
            
            disp(['Not found: ',this_strain,'_WT']);
            
        end

    end
    
    params.colors=params_colors;
    params.strains=params_strains;
    params.species=params_species;
    params.cs=params_cs;
    params.Vs=params_Vs;
    params.Ks=ones(1,2*numStrains); %*1e-9; <--???
    
    %params.cs=[bp_cs, bp_cs];
    %params.Vs=[(params.vcost).*bp_Vs, bp_Vs];
    %params.Ks=[bp_Ks,bp_Ks];
    
    

    params.S0=S0;
    params.T=T;
    params.d=d;

    params.epsilon=epsilon;
    params.numStrains=numStrains;
