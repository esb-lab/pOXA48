# pNUK73: A Metropolis-Hastings MCMC implementation to fit a Monod model to a bacterial growth curve
# Version: 1.0.1
# Date: 2014-07-31
#
# Author Rafael Pena-Miller <rpm@ccg.unam.mx> based in an MCMC implementation by Ben S. Cooper
# Copyright 2014. All rights reserved.
#
# This script receives parameters from command line and executes MH_MCMC.R
# 
# Example: 
# $ Rscript runMCMC.R dataFile="../data/OD600_B.csv" w.prior="uniform" N=1e6
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
# Receives and parses arguments from command line
args <- commandArgs(trailingOnly = TRUE)
if(length(args)>0){
  for (i in 1:length(args) ) {
    arg<-unlist(strsplit(args[i], "="))
    if(arg[1]=="dataFile"){
      dataFile<-arg[2]
    }else if(arg[1]=="runDir"){
      runDir<-arg[2]
    }else if(arg[1]=="w.prior"){
      w.prior<-arg[2]
    }else if(arg[1]=="n.clones"){
      n.clones<-as.numeric(arg[2])
    }else if(arg[1]=="N"){
      N<-as.numeric(arg[2])
    }
  }
}


if(!exists("dataFile"))  stop("dataFile not defined")
if(!exists("w.prior"))  w.prior<-"uniform"
if(!exists("n.clones"))  n.clones<-1
if(!exists("N"))  N<-1e5
if(!exists("runDir")) runDir=paste("../runs/_",format(Sys.time(), "%Y%m%d_%H%M%S"),"_",w.prior,sep="")

#Create output directory if not exists
if (!file.exists(runDir)){
  dir.create(runDir,recursive=TRUE)
  sink(paste(runDir,'/summary.txt', sep = ""), append=TRUE, split=TRUE) 
}

print(noquote(paste("Output Directory: ",runDir)))
print(noquote(paste("Input Data file: ",dataFile)))
print(noquote(paste("Number of iterations: ",N)))
print(noquote(paste("Prior: ",w.prior)))
print(noquote(paste("Number of clones: ",n.clones)))

############################################
# Load data 
ODdata<-read.csv(dataFile)

############################################
# Clone data 
ODdata <- lapply(as.data.frame(ODdata), function(z) rep(z, n.clones))
ODdata <- as.data.frame(ODdata)
ODdata<-ODdata[order(ODdata$t),]

############################################
# Runs MCMC script
start.time <- Sys.time()
print(noquote("*******************************"))

source("MH_MCMC.R")

end.time <- Sys.time()
time.taken <- end.time- start.time
print(noquote("MCMC Done...")) # Execution time
print(noquote(time.taken)) # Execution time

############################################
# Save variables into a file
save.image(file=paste(runDir,"/MCMC.RData", sep=""))
print(noquote(paste("image saved in ",paste(runDir,"/MCMC.RData", sep="")))) # Execution time

############################################
# Produce MCMC diagnostic plots
#source("plotMCMC.R")

print(noquote(" Export done."))

#print(posterior.samples$rho)
#print(posterior.samples$muK)

write.csv(posterior.samples$rho,paste(runDir,"posterior_samples_rho.csv", sep = ""), row.names = FALSE)
write.csv(posterior.samples$muK,paste(runDir,"posterior_samples_muK.csv", sep = ""), row.names = FALSE)

