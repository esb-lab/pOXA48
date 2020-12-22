

function params=randomParameters(numStrains, T, d, S0, seg_rate, conj_rate, epsilon, mean_muKs, std_muKs, mean_rhos, std_rhos)

    params.seg_rate=seg_rate;
    params.conj_rate=conj_rate;


    params.S0=S0;
    params.T=T;
    params.d=d;

    params.epsilon=epsilon;
    params.numStrains=numStrains;
    
    %Growth parameters
    cs_TC=mean_rhos(1).*abs(1+normrnd(0,std_rhos(1)/mean_rhos(1), [1, numStrains]));
    Vs_TC=mean_muKs(1).*abs(1+normrnd(0,std_muKs(1)/mean_muKs(1), [1, numStrains]));
    
    cs_WT=mean_rhos(2).*abs(1+normrnd(0,std_rhos(2)/mean_rhos(2), [1, numStrains]));
    Vs_WT=mean_muKs(2).*abs(1+normrnd(0,std_muKs(2)/mean_muKs(2), [1, numStrains]));
    
    params.cs=[cs_TC, cs_WT];
    params.Vs=[Vs_TC, Vs_WT];
    params.Ks=ones(1,length(params.Vs));

    params.strains={};
    params.species={};
    params.plasmids={};
    for i=1:numStrains*2
        params.strains{i}=['strain',num2str(i)];
        params.species{i}='sim';
        if i<length(cs_TC)
            params.plasmids{i}='TC';
        else
            params.plasmids{i}='WT';
        end
    end