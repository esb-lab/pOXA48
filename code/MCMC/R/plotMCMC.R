# pNUK73: A Metropolis-Hastings MCMC implementation to fit a Monod model to a bacterial growth curve
# Version: 1.0.1
# Date: 2014-07-31
#
# Author Rafael Pena-Miller <rpm@ccg.unam.mx> based in an MCMC implementation by Ben S. Cooper
# Copyright 2014. All rights reserved.
#
# This script plots MCMC diagnostic plots based on the output of MH_MCMC.R
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

library(Hmisc)
library(extrafont)
library(car)



# Create directory to save plots
if (!file.exists(paste(runDir,"/figures/",sep=""))){
  dir.create(paste(runDir,"/figures/",sep=""), recursive=TRUE)
}

##############################################
# This function computes the variances (and means) of the data for the 
# intervals from 1 to incr, 1 to 2*incr, 1 to 3*incr, ..., till finally 
# 1 to length(x) is calculated. The lengths of the corresponding intervals 
# is given in the index. 
# Colin J. Greene, July 11, 1995     modified by JM Hoenig June 2007 to run in R
updatevar <- function( x, incr=length(x) ){
  lx <- length(x)
  if( lx < incr ){
    stop("Increment must be less than or equal to the length of data") }
  else if( lx == incr ){
    updatevals <- vector("numeric", length=0) }
  else if( lx < 2*incr ){
    updatevals <- c( lx ) }
  else if( (lx %% incr) == 0 ){
    updatevals <- seq(2*incr,lx,by=incr) }
  else {
    updatevals <- c( seq(2*incr,lx,by=incr), lx) }
  j <- p1 <- p2 <- p3 <- p4 <- nmean <- nvar <- 1
  nmean[j] <- mean(x[1:incr])
  nvar[j]  <- var(x[1:incr])
  nval     <- incr
  for(k in updatevals){
    nmean[j+1] <- mean(x[1:k])
    p1 <- (nval-1)*nvar[j]
    p2 <- sum( (x[(nval+1):k]-nmean[j])^2 )
    p3 <- 2*(nmean[j]-nmean[j+1])*sum( (x[(nval+1):k]-nmean[j]) )
    p4 <- (k)*(nmean[j]-nmean[j+1])^2
    nvar[j+1] <- (p1 + p2 + p3 + p4)/(k-1)
    j <- j + 1
    nval <- k
  }
  return(list( var=nvar,mean=nmean,index=c(incr,updatevals) ))
  invisible()
}

################################################
#Plots running means and variances
Nx<-500
inix<-list(rho=5e8, muK=1e-10, disp=1e8) 

#Plots running mean and variance for rho
xrho<-posterior.samples$rho[1:Nx]
meanxrho<-updatevar(xrho,2)
vardownrho<-meanxrho$mean - sqrt(meanxrho$var)
varuprho<-meanxrho$mean + sqrt(meanxrho$var)

pdf(paste(runDir,'/figures/running_mean_rho.pdf',sep = ""))
par(mar=c(5.1,5.1,4.1,2.1))
plot(c(posterior.samples$rho[1], meanxrho$mean),type="l",main=expression(rho), ylim=c(0.975*min(vardownrho),1.01*max(varuprho)), xlab="Iteration", ylab="Posterior mean", col="white",cex.lab=2,cex.main=2, cex.axis = 1.5)
polygon(c(2:(length(varuprho)+1), (length(vardownrho)+1): 2), c(varuprho,rev(vardownrho)),col="gray96",border=NA)
lines(c(inix$rho, varuprho),col="gray")
lines(c(inix$rho, vardownrho),col="gray")
lines(c(inix$rho, meanxrho$mean),type="l",main=expression("rho"),lw=3)
dev.off()

######
#Plots running mean and variance for rho
xmuK<-posterior.samples$muK[1:Nx]
meanxmuK<-updatevar(xmuK,2)
vardownmuK<-meanxmuK$mean - sqrt(meanxmuK$var)
varupmuK<-meanxmuK$mean + sqrt(meanxmuK$var)

pdf(paste(runDir,'/figures/running_mean_muK.pdf',sep = ""))
par(mar=c(5.1,5.1,4.1,2.1))
plot(c(posterior.samples$muK[1], meanxmuK$mean),type="l",main=expression(bar(mu)/K), ylim=c(0.975*min(vardownmuK),1.01*max(varupmuK)), xlab="Iteration", ylab="Posterior mean", col="white",cex.lab=2,cex.main=2, cex.axis = 1.5)
polygon(c(2:(length(varupmuK)+1), (length(vardownmuK)+1): 2), c(varupmuK,rev(vardownmuK)),col="gray96",border=NA)
lines(c(inix$muK, varupmuK),col="gray")
lines(c(inix$muK, vardownmuK),col="gray")
lines(c(inix$muK, meanxmuK$mean),type="l",main=expression("muK"),lw=3)
dev.off()

######
# Plots 2D distribution 
png(paste(runDir,'/figures/muK_rho.png', sep = ""),width = 3150, height = 3150, units = "px",pointsize = 12, bg = "white", res=300, family = "sans")
par(mar=c(5.1, 5, 5.1, 1))
plot(posterior.samples$rho,posterior.samples$muK,xlab=expression(rho),ylab=expression(bar(mu)/K),type='p',pch=".",cex.lab=2,cex.main=2, cex.axis = 1.25, family = "sans")
dev.off()

######
# Plots histograms
pdf(paste(runDir,'/figures/histogram_rho.pdf',sep = ""))
hist(posterior.samples$rho, breaks=25, col="white")
dev.off()

pdf(paste(runDir,'/figures/histogram_muK.pdf',sep = ""))
hist(posterior.samples$muK, breaks=25, col="white")
dev.off()

######
# Plots chains
png(paste(runDir,'/figures/chains.png', sep = ""),width = 2100, height = 2100, units = "px",pointsize = 12, bg = "white", res=300, family = "sans")
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(1,1), respect = FALSE)
par(mar=c(1, 5, 5.1, 1))

plot(posterior.samples$muK,type='l',xlim=c(0, length(posterior.samples$muK)), ylab = expression(bar(mu)/K),xaxt="n",cex.lab=2,cex.main=2, cex.axis = 1.2, family = "sans")
axis(side = 1, tck = -.05, labels = NA)
#text(90,u[4],labels = expression(rho),col = "black",pos = 1,cex=2)

par(mar=c(5.1, 5, 1, 1))
plot(posterior.samples$rho,type='l',xlim=c(0, length(posterior.samples$muK)), ylab = expression(rho),cex.lab=2,cex.main=2, cex.axis = 1.2, family = "sans")
dev.off()

######
# Plots chain autocorrelation
pdf(paste(runDir,'/figures/autocorrelation.pdf', sep = ""),width = 7, height = 7,pointsize = 12, bg = "white", useDingbats=FALSE, family="sans")
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(1,1), respect = FALSE)
par(mar=c(1, 5, 5.1, 1))

r <- acf(posterior.samples$muK,plot=F, lag.max=100)
plot(r$lag[2:length(r$acf)], r$acf[2:length(r$acf)], type = "h", lwd = 2, col = "black", ylab = "",xaxt="n",cex.lab=2,cex.main=2, cex.axis = 1.25, family = "sans")
axis(side = 1, tck = -.05, labels = NA)
ci <- .95 
clim <- qnorm( (1+ci) / 2 ) / sqrt(r$n.used)
abline(h = c(-1,1) * clim, lty = 2, col = "black", lwd = .5)

par(mar=c(5.1, 5, 1, 1))
r <- acf(posterior.samples$rho,plot=F, lag.max=100)
plot(r$lag[2:length(r$acf)], r$acf[2:length(r$acf)], type = "h", lwd = 2, col = "black",xlab = "Lag", ylab = "",cex.lab=2,cex.main=2, cex.axis = 1.25, family = "sans")
ci <- .95 
clim <- qnorm( (1+ci) / 2 ) / sqrt(r$n.used)
abline(h = c(-1,1) * clim, lty = 2, col = "black", lwd = .5)
#axis(side = 2, tck = -.05, at = c(-round(clim,3),0,round(clim,3)),cex.axis=1.25, family = "sans")
u <- par("usr")
#text(90,u[4],labels = expression(bar(mu)/K),col = "black",pos = 1,cex=2)
mtext("Autocorrelation", line=3, side=2, cex=2, at=u[4], family = "sans")
dev.off()

######
# Plots 6 randomly selected runs
pdf(paste(runDir,'/figures/selected_fits.pdf',sep = ""),width = 7, height = 7,pointsize = 12, bg = "white", useDingbats=FALSE, family="sans")
par(mfrow=c(3,2),mar=c(5,5,1,1))
selected.runs<-  c(1:n.samples)[order(runif(n.samples))][1:6]
for(i in selected.runs){
  parameters<-c(rho=as.numeric(posterior.samples[i,]$rho), mu=this_mu, K=this_mu/as.numeric(posterior.samples[i,]$muK))
  selected_out<-ode(times=times, y=yini, func="derivs", parms=parameters, jacfun="jac", dllname="ODfit", initfunc="initmod",nout=1, outnames="OD")
  errbar(ODdata[,1], OBSERVED.OD, OBSERVED.OD-OBSERVED.stOD, OBSERVED.OD+OBSERVED.stOD, xlab='Time (hours)', ylab='Bacterial Density',cex.lab=2,cex.main=2, cex.axis = 1.25, family = "sans")
  lines(selected_out[,1],(selected_out[,3]),col='red', lwd=2)
  axlim=par("usr")
  text(.4*axlim[2],axlim[3]+1.5*(axlim[4]-axlim[3])/10,sprintf('rho=%.4e',posterior.samples[i,]$rho),pos=4,cex=1.5)
  text(.4*axlim[2],axlim[3]+3*(axlim[4]-axlim[3])/10,sprintf('mu/K=%.4e',posterior.samples[i,]$muK),pos=4,cex=1.5)
  if(i==selected.runs[1]){
    legend(0,1.2e9, c("Data","Model"), lty=c(1,1), lwd=c(2.5,2.5),col=c("black","red"),cex=1.25) 
  }
}
dev.off()

######
# Plots of predicted v observed 
P<-100  # First sample P runs for plotting randomly selected without replacement
runs.to.plot<-c(1:n.samples)[order(runif(n.samples))][1:P]
predicted<-NULL
residual<-NULL
for(i in runs.to.plot){
  parameters<-c(rho=as.numeric(posterior.samples[i,]$rho), mu=this_mu, K=this_mu/as.numeric(posterior.samples[i,]$muK))
  sampled_out<-ode(times=times, y=yini, func="derivs", parms=parameters, jacfun="jac", dllname="ODfit", initfunc="initmod",nout=1, outnames="OD")
  xdat=seq(from=0,to=sampled_out[length(sampled_out[,3]),1],length.out=length(OBSERVED.OD))
  Bpred=pchip(sampled_out[,1], sampled_out[,3], xdat)
  predicted<-rbind(predicted,Bpred)
}

mean_predicted<-apply(predicted, 2, mean)
residual<-mean_predicted-OBSERVED.OD

######
# Plots predicted vs observed (qq-plot)
pdf(paste(runDir,'/figures/qq_predicted_observed.pdf', sep = ""),pointsize = 12, bg = "white")
qqplot(mean_predicted,OBSERVED.OD,xlab="Theoretical quantiles",ylab="Sample quantiles",type="p", main="Normal Q-Q plot", cex = 1.5, lwd=2)
abline(0,1, col = "black", lty=2)
dev.off()

######
# Plots OD-dependant residuals
pdf(paste(runDir,'/figures/residuals_OD.pdf', sep = ""),pointsize = 12, bg = "white",width = 7, height = 7,useDingbats=FALSE, family="sans")
par(mar=c(5.1, 5, 5.1, 1))
plot(mean_predicted,residual,xlab="Predicted bacterial density",ylab="Residual",type="p", col="black",bg="black",lwd=3,cex.lab=2,cex.main=2, cex.axis = 1.25, family = "sans")
abline(h = 0, col = "gray60", lty=2)
dev.off()


######
# Plots time-dependant residuals 
pdf(paste(runDir,'/figures/residuals_time.pdf', sep = ""),pointsize = 12, bg = "white")
plot(xdat,residual,xlab="Hours",ylab="Residual",type='l',pch=".", lwd=2)
abline(h = 0, col = "gray60", lty=2)
dev.off()

######
# Plots residual distribution
pdf(paste(runDir,'/figures/residuals_distribution.pdf',sep = ""))
hist(residual, breaks=15, col="white")
dev.off()

######
# Plots residuals (qqplot)
pdf(paste(runDir,'/figures/residuals_qq.pdf', sep = ""),pointsize = 12, bg = "white")
qqPlot(residual)
dev.off()


print(noquote(" Plots done."))


