function [muK_R, rho_R]=computeGelmanRubin(dataDir, M)

    rho_theta_m=zeros(1,M);
    rho_sigma2_m=zeros(1,M);
    
    muK_theta_m=zeros(1,M);
    muK_sigma2_m=zeros(1,M);
    
    N_m=zeros(1,M);

    for num_chain=1:M

        if exist([dataDir, 'posterior_samples_rho_num_chain',num2str(num_chain),'.csv'],'file')

            this_posterior_samples_rho = table2array(readtable([dataDir, 'posterior_samples_rho_num_chain',num2str(num_chain),'.csv']));
            this_posterior_samples_muK = table2array(readtable([dataDir, 'posterior_samples_muK_num_chain',num2str(num_chain),'.csv']));

            rho_theta_m(num_chain)=mean(this_posterior_samples_rho);
            rho_sigma2_m(num_chain)=sum((rho_theta_m(num_chain)-this_posterior_samples_rho).^2)/(length(this_posterior_samples_rho)-1);
            
            muK_theta_m(num_chain)=mean(this_posterior_samples_muK);
            muK_sigma2_m(num_chain)=sum((muK_theta_m(num_chain)-this_posterior_samples_muK).^2)/(length(this_posterior_samples_rho)-1);

            
            N_m(num_chain)=length(this_posterior_samples_rho);


        end
    end



    %Calculate theta, the mean of all chains
    rho_theta=mean(rho_theta_m);
    muK_theta=mean(muK_theta_m);

    %Compute how the individual means vary around the joint mean
    rho_B=0;
    muK_B=0;
    for m=1:M
       rho_B=rho_B+(N_m(m)/(M-1))*(rho_theta_m(m)-rho_theta)^2;
       muK_B=muK_B+(N_m(m)/(M-1))*(muK_theta_m(m)-muK_theta)^2;
    end

    %Compute the averaged variances of the chains
    rho_W=(1/M)*sum(rho_sigma2_m);
    muK_W=(1/M)*sum(rho_sigma2_m);
    
    N=mean(N_m);
    
    rho_V=((N-1)/N)*rho_W+ ((M+1)/(M*N))*rho_B;
    muK_V=((N-1)/N)*muK_W+ ((M+1)/(M*N))*muK_B;

    rho_R=sqrt(rho_V/rho_W);
    muK_R=sqrt(muK_V/muK_W);
    
    