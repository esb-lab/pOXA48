function fout = fMany(~, y, params)

    S=y(1);
    numStrains=(length(y)-1)/2;

    B_TC=y(2:numStrains+1);
    B_WT=y(numStrains+2:end);

    uStot=0;
    dB=zeros(1, numStrains);

    %For plasmid-bearing (TC)
    for i=1:numStrains
        uSi=uS(S, params.Vs(i), params.Ks(i));
        uStot=uStot+uSi*B_TC(i);  %Resource consumed

        %Plasmid transfer B_TC_j -> B_WT_i
        dBconj=0;
        for j=1:numStrains
                dBconj=params.conj_rate*B_TC(i)*B_WT(j);
        end
        dB(i)=params.cs(i)*uSi*B_TC(i)*(1 - params.seg_rate) + dBconj  - params.d*B_TC(i);
    end

    %For plasmid-free (WT)
    for i=1:numStrains
        uSi=uS(S, params.Vs(numStrains+i), params.Ks(numStrains+i));
        uStot=uStot+uSi*B_WT(i);  %Resource consumed

        %Plasmid transfer B_TC_j -> B_WT_i
        dBconj=0;
        for j=1:numStrains
                dBconj=params.conj_rate*B_TC(j)*B_WT(i);
        end
        dB(numStrains+i)=params.cs(numStrains+i)*uSi*B_WT(i) - dBconj + params.cs(i)*uSi*B_TC(i)*(params.seg_rate) - params.d*B_WT(i);
    end

    dS=-uStot + params.d*(params.S0-S);
    fout = params.T*([dS dB]');
end


function ret = uS(S,V,K)
    ret=(S*V)/(K+S);
end
