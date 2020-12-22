function gm=getGMModel(MCMC_muKs, MCMC_rhos, max_muKs, max_rhos)

if nargin <3
    max_muKs=10e-10;
    max_rhos=12e8;
end


nMCMC_muKs=MCMC_muKs./max_muKs;
nMCMC_rhos=MCMC_rhos./max_rhos;

iok=intersect(find(MCMC_muKs<max_muKs) ,find(MCMC_rhos<max_rhos)) ;
X = [nMCMC_muKs(iok)' nMCMC_rhos(iok)'];
GMModel = fitgmdist(X,1);
nfitMu=GMModel.mu;
nfitSigma=GMModel.Sigma;
%fitMu=nfitMu.*[max_muKs max_rhos];
%fitSigma=nfitSigma.*[max_muKs max_rhos];
%gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1)); g = gca;

gm = gmdistribution(nfitMu,nfitSigma);
