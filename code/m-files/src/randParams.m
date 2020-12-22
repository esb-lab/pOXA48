

function params=randParams(numStrains, T, d, S0, meanCost, cstar, Vstar, Kstar, seg_rate, conj_rate, epsilon, sigma_cost, sigma_growth)


    %Growth parameters
    bp_cs=cstar.*abs(1+normrnd(0,sigma_growth(1), [1, numStrains]));
    bp_Vs=Vstar.*abs(1+normrnd(0,sigma_growth(2), [1, numStrains]));
    bp_Ks=Kstar.*abs(1+normrnd(0,sigma_growth(3), [1, numStrains]));

    %params.cost=cost;
    params.vcost=meanCost+normrnd(0,sigma_cost,[1, numStrains]);

    params.seg_rate=seg_rate;
    params.conj_rate=conj_rate;

    params.cs=[(params.vcost).*bp_cs, bp_cs];
    params.Vs=[(params.vcost).*bp_Vs, bp_Vs];

    params.Ks=[bp_Ks,bp_Ks];

    params.S0=S0;
    params.T=T;
    params.d=d;

    params.epsilon=epsilon;
    params.numStrains=numStrains;
    
    
    params.strains={};
    params.species={};
    params.plasmids={};
    for i=1:numStrains
        params.strains{i}=['strain',num2str(i)];
        params.strains{i+numStrains}=['strain',num2str(i)];
        
        params.species{i}='sim';
        params.species{i+numStrains}='sim';
        
        params.plasmids{i}='TC';
        params.plasmids{i+numStrains}='WT';
    end

