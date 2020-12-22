

function params=gmParameters(gm_sample,  T, d, S0, seg_rate, conj_rate, epsilon)

    numStrains=length(gm_sample)/2;
    
    if numStrains<1e2
        cmap0 = distinguishable_colors(numStrains);
    else
        cmap0=ones(numStrains,3)*.5;
    end
    cmap=[cmap0; cmap0];

    %Growth parameters

    params.seg_rate=seg_rate;
    params.conj_rate=conj_rate;

    params.S0=S0;
    params.T=T;
    params.d=d;

    params.epsilon=epsilon;
    params.numStrains=numStrains;
    
    params_species={};
    params_plasmids={};
    params_strains={};
    params_colors=zeros(2*numStrains,3);
    params_cs=zeros(1,2*numStrains);  %rhos
    params_Vs=zeros(1,2*numStrains);  %Vmax
    params_Ks=zeros(1,2*numStrains);  %Km
    for icomm=1:numStrains
       
            %Plasmid-bearing (TC)
            params_cs(icomm)=gm_sample(icomm,2);
            params_Vs(icomm)=gm_sample(icomm,1);
            params_Ks(icomm)=1; %1e-9;
            params_species{icomm}='ODE';
            params_strains{icomm}=['strain',num2str(icomm)];
            params_plasmids{icomm}='TC';
            params_colors(icomm,:)=cmap(icomm,:); 
            
            %Plasmid-bearing (TC)
            params_cs(icomm+numStrains)=gm_sample(icomm+numStrains,2);
            params_Vs(icomm+numStrains)=gm_sample(icomm+numStrains,1);
            params_Ks(icomm+numStrains)=1; %1e-9;
            params_species{icomm+numStrains}='ODE';
            params_strains{icomm+numStrains}=['strain',num2str(icomm)];
            params_plasmids{icomm+numStrains}='WT';
            params_colors(icomm+numStrains,:)=cmap(icomm,:); 
            
    end
    params.cs=params_cs;
    params.Vs=params_Vs;
    params.Ks=params_Ks;
    params.species=params_species;
    params.plasmids=params_plasmids;
    params.strains=params_strains;
    params.colors=params_colors;
    