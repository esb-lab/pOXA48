function sub_params=subsetParameters(comm_params, N)
    
    p = randperm(comm_params.numStrains,N);
    
    sub_params.index=p;
    
    sub_params.vcost=comm_params.vcost(p);

    sub_params.seg_rate=comm_params.seg_rate;
    sub_params.conj_rate=comm_params.conj_rate;

    sub_params.cs=comm_params.cs([p p+comm_params.numStrains]);
    sub_params.Vs=comm_params.Vs([p p+comm_params.numStrains]);

    sub_params.Ks=comm_params.Ks([p p+comm_params.numStrains]);

    sub_params.S0=comm_params.S0;
    sub_params.T=comm_params.T;
    sub_params.d=comm_params.d;

    sub_params.epsilon=comm_params.epsilon;
    sub_params.numStrains=N;
