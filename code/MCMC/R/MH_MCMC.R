# pNUK73: A Metropolis-Hastings MCMC implementation to fit a Monod model to a bacterial growth curve
# Version: 1.0.1
# Date: 2014-07-31
#
# Author Rafael Pena-Miller <rpm@ccg.unam.mx> based in an MCMC implementation by Ben S. Cooper
# Copyright 2014. All rights reserved.
#
# pNUK73 provides a simple MCMC method to estimate growth kinetic parameters (mu, Ks, rho) 
# that fit growth curve measured as bacterial optical density.
#
# In particular we want to make inference about the following parameters: 
# bar(mu), K, rho: Basic growth kinetic parameters of a single bacterial type growing under a single-limiting resource
# bar(mu) - Maximum uptake rate
# K - Half saturation constant
# rho - resource conversion rate
# Note: bar(mu) and K are not identifiable, but their ratio is 
# Therefore we would fit a two-parameter Monod model: resource conversion rate (rho) and specific affinity (bar(mu)/K)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
############################################

library(deSolve)
library(pracma)

############################################
# MCMC default parameters

if(!exists("ODdata"))  stop("ODdata not defined")
if(!exists("N"))  N<-1e6
if(!exists("thin"))  thin<-1e2
if(!exists("fracBurnIn"))  fracBurnIn=.2 #in percentage of runs
if(!exists("target.acceptance"))  target.acceptance<-.23
if(!exists("target.tolerance"))  target.tolerance<-.001
if(!exists("w.prior"))  w.prior<-"uniform"
if(!exists("n.clones"))  n.clones<-1


############################################
# Function that converts OD600 into CFUs 
# Based on an experiment by ASM
CFUs<-function(ODs){
  b1<-1002614811.71;
  b2<-12526559.79;
  CFUs<-b1*ODs+b2;
  return(CFUs)
}

OBSERVED.OD<-CFUs(ODdata[,2])
OBSERVED.stOD<-CFUs(ODdata[,3])
OD0=CFUs(ODdata[1,2])  #Initial condition used in the simulations from data
this_mu<-1e-9  #Fix one of the parameters so the other will be given by the determined ratio

############################################
# # Model: One bacterial type, single limiting resource
# Simulate model using deSolve library (from a compiled version)
system("R CMD SHLIB ../src/ODfit.c")
dyn.load("../src/ODfit.so")

# ODE Params
times <- seq(from=0,to=tail(ODdata[,1],n=1),by=.5)
yini <- c(R=1, B=OD0)  #Initial condition from data

############################################
ODfit.mh<-function(parameters, yini, i=0){
  #The version is designed for use with Metropolis Hastings algorithm
  #yini represents the initial OD
  with(as.list(c(parameters, yini)), {
     
    #define K from mu/K ratio and a fixed value of mu
    muK=as.numeric(parameters['muK'])
    this_K=this_mu/muK
    parameters1<-c(rho=as.numeric(parameters['rho']), mu=this_mu, K=this_K)
    
    #Simulate model using deSolve compiled model
    out<-ode(times=times, y=yini, func="derivs", parms=parameters1, jacfun="jac", dllname="ODfit", initfunc="initmod",nout=1, outnames="OD")
    
    #Both the observed and the predicted vectors have to be the same length
    xdat=seq(from=0,to=out[length(out[,3]),1],length.out=length(OBSERVED.OD))
    Bpred=pchip(out[,1], out[,3], xdat)
    
    # now construct a likelihood  
    # below OBSERVED.OD is the actual data (specified as global variable)
    LLterms<-dnorm(OBSERVED.OD, mean=Bpred, sd=as.numeric(parameters['disp']), log=TRUE)
    LL<-n.clones*sum(LLterms)  #Multiplied by the number of clones
  
    return(list(LL=LL))
  })   
}

############################################
###### Metropolis-Hastings sampling
############################################

# 1. Define Prior distributions
if(w.prior=="uniform"){
  # Uniform prior distribution
  muK.prior.density<-function(x){ dunif(x, min=0, max=1e-8)}  
  rho.prior.density<-function(x){ dunif(x, min=0, max=1e11)}  
  disp.prior.density<-function(x){ dunif(x, min=0.1, max=1e13)}
  
}else if(w.prior=="lognormal"){
  # Lognormal prior distribution
  muK.prior.density<-function(x){ dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)}  
  rho.prior.density<-function(x){ dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)}  
  disp.prior.density<-function(x){ dunif(x, min=0.1, max=1e13)}

}else if(w.prior=="gamma"){
  # Gamma prior distribution
  muK.prior.density<-function(x){ dgamma(x, shape=0.001,rate=0.001)}
  rho.prior.density<-function(x){ dunif(x, min=0, max=1e11)}  
  disp.prior.density<-function(x){ dunif(x, min=0.1, max=1e13)}

}else if(w.prior=="beta"){
  # Beta prior distribution
  muK.prior.density<-function(x){ dbeta(x, 1,1)}
  rho.prior.density<-function(x){ dunif(x, min=0, max=1e11)}  
  disp.prior.density<-function(x){ dunif(x, min=0.1, max=1e13)}
}

######
# 2. Define functions to propose new values for parameters

# Initial variance of proposals 
std<-list(rho=5e7, muK=5e-11, disp=5e7)



# Adjust (adaptively) variance of the proposals to tune the algorithm (ideally we want to accept about a fifth of proposals for each variable)
# Note that these functions may propose values outside possible ranges, but since prior distributions
# will have zero probability at these values they have no chance of being accepted
generate.proposal.muK<-function(currentmuK, currentsd){return(currentmuK + rnorm(1,mean=0,sd=currentsd))}
generate.proposal.rho<-function(currentrho, currentsd){return(currentrho + rnorm(1,mean=0,sd=currentsd))}
generate.proposal.disp<-function(currentdisp, currentsd){return(currentdisp + rnorm(1,0,sd=currentsd))}

######
#run first time
parameters<-list(rho=5e8, muK=5e-10, disp=5e8) 
proposed.parameters<-parameters
current.loglikelihood<-ODfit.mh(parameters, yini, 0)



#store results from the posterior sample with one column for each parameter & one for LL & one for res
posterior.samples<-as.data.frame(matrix(c(as.numeric(parameters),current.loglikelihood$LL),nrow=1))
colnames(posterior.samples)<-c(names(as.vector(parameters)),"LL")

######

# it will run without saving values until every parameter is inside the tolerated interval
hit<-list(rho=0, muK=0, disp=0) 

#variables that store the number of accepted moves for each variable
rho.moves.accepted<-0
muK.moves.accepted<-0
disp.moves.accepted<-0

for(i in 1:N){
  
  ######
  #  updates to muK
  proposed.parameters$muK<-generate.proposal.muK(parameters$muK, std$muK) # proposed update
  current.prior.density<-muK.prior.density(parameters$muK)
  proposal.prior.density<-muK.prior.density(proposed.parameters$muK)
    if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
      proposal.loglikelihood<-ODfit.mh(proposed.parameters, yini, i)
      acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
      if(runif(1)<acceptance.prob){  
        parameters <- proposed.parameters  
        current.loglikelihood<- proposal.loglikelihood
        muK.moves.accepted<-muK.moves.accepted+1
      } else proposed.parameters<-parameters  # otherwise current parameters stay the same
    } else  proposed.parameters<-parameters
  
  ######
  # updates to rho
  proposed.parameters$rho<-generate.proposal.rho(parameters$rho, std$rho) # proposed update 
  current.prior.density<-rho.prior.density(parameters$rho)
  proposal.prior.density<-rho.prior.density(proposed.parameters$rho)
    if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
      proposal.loglikelihood<-ODfit.mh(proposed.parameters, yini,  i) 
      acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL) * proposal.prior.density/current.prior.density)
      
      if(runif(1)<acceptance.prob){ 
        parameters <- proposed.parameters  
        current.loglikelihood<- proposal.loglikelihood
        rho.moves.accepted<-rho.moves.accepted+1
      }  else proposed.parameters<-parameters # otherwise current parameters stay the same
    } else  proposed.parameters<-parameters
  
  #  updates to disp (dispersion parameter for negative binomial likeihood [if it is NA use Poisson likelihood insted])
  if(!is.na(parameters$disp)){
    proposed.parameters$disp<-generate.proposal.disp(parameters$disp, std$disp) # proposed update 
    current.prior.density<-disp.prior.density(parameters$disp)
    proposal.prior.density<-disp.prior.density(proposed.parameters$disp)
      if(proposal.prior.density>0) { # no point running the model if proposal is outside feasbile range as acceptance.prob will be zero
        proposal.loglikelihood<-ODfit.mh(proposed.parameters, yini,  i) 
        acceptance.prob<-min(1, exp(proposal.loglikelihood$LL- current.loglikelihood$LL)* proposal.prior.density/current.prior.density)
        
        if(runif(1)<acceptance.prob){ 
          parameters <- proposed.parameters  
          current.loglikelihood<- proposal.loglikelihood
          disp.moves.accepted<-disp.moves.accepted+1
        }  else proposed.parameters<-parameters # otherwise current parameters stay the same
      } else  proposed.parameters<-parameters
    }else proposed.parameters<-parameters
  
  ######
  
  if(i%%(N/10)==0){  #Display progress...
    print(noquote(paste("...",100*i/N,"%",sep="")))
  }
  
  # now store results of this iteration if i is a multiple of the parameter "thin"
  if(i%% thin ==0) {
    if(hit$rho && hit$muK && hit$disp){ 
      #if the acceptance rate of each parameter is inside the tolerated interval
      new.posterior.sample<-as.data.frame(matrix(c(as.numeric(parameters),current.loglikelihood$LL),nrow=1))
      colnames(new.posterior.sample)<-c(names(as.vector(parameters)),"LL")
      posterior.samples<-rbind(posterior.samples, new.posterior.sample)
      
    }else{
      #correct variance of proposals in order to accept about a fifth of each variables
      # (and do not save the samples)
      #For muK
        this.acceptance<-muK.moves.accepted/i
        if(this.acceptance<target.acceptance-target.tolerance){
          std$muK<-std$muK-std$muK/10
        }else if(this.acceptance>target.acceptance+target.tolerance){
          std$muK<-std$muK+std$muK/10
        }else{
          hit$muK=T
        }
      
      #For rho
        this.acceptance<-rho.moves.accepted/i
        if(this.acceptance<target.acceptance-target.tolerance){
          std$rho<-std$rho-std$rho/10
        }else if(this.acceptance>target.acceptance+target.tolerance){
          std$rho<-std$rho+std$rho/10
        }else{
          hit$rho=T
        }
      
      #For disp
        this.acceptance<-disp.moves.accepted/i
        if(this.acceptance<target.acceptance-target.tolerance){
          std$disp<-std$disp-std$disp/10
        }else if(this.acceptance>target.acceptance+target.tolerance){
          std$disp<-std$disp+std$disp/10
        }else{
          hit$disp=T
        } 
      
    }
    
  }
} #  end of main MCMC loops

# define burnIn period based on the number of samples and the fraction to chuck out
burnIn=floor(length(posterior.samples$muK)*fracBurnIn)

# now chuck out first burnIn samples (i.e. burnIn*thin iterations) as a burnin 
posterior.samples<-posterior.samples[burnIn: dim(posterior.samples)[1],]
n.samples<-length(posterior.samples[,1])

############################################
#  Summzarize results
print(noquote(" "))
print(noquote(sprintf("N=%0.1e", N)))
print(noquote(sprintf("fracBurnIn=%0.2f", fracBurnIn)))
print(noquote(sprintf("thin=%0.0f", thin)))
if(n.clones>1){ #If cloning data
  print(noquote(sprintf("K=%0.0f", n.clones)))
}

# calc stats  - mean and 95% CIs for key quantities
# 1. rho
mean_rho=mean(posterior.samples$rho)
var_rho=var(posterior.samples$rho)
std_rho=std(posterior.samples$rho)
print(noquote(" "))
print(noquote(sprintf("** mean(rho)=%e", mean_rho)))
print(noquote(sprintf("** var(rho)=%e", var_rho)))
print(noquote(sprintf("** std(rho)=%e", std_rho)))
print(noquote(quantile( posterior.samples$rho, c(.025,.5,.975))))

# 2. muK
mean_muK=mean(posterior.samples$muK)
var_muK=var(posterior.samples$muK)
std_muK=std(posterior.samples$muK)
print(noquote(" "))
print(noquote(sprintf("** mean(muK)=%e", mean_muK)))
print(noquote(sprintf("** var(muK)=%e", var_muK)))
print(noquote(sprintf("** std(muK)=%e", std_muK)))
print(noquote(quantile(posterior.samples$muK, c(.025,.5,.975))))
print(noquote(" "))






