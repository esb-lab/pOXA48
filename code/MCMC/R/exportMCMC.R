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

print(paste(">",runDir))

print(noquote(" Export done."))

#print(posterior.samples$rho)
#print(posterior.samples$muK)

write.csv(posterior.samples$rho,paste(runDir,"posterior_samples_rho.csv", sep = ""), row.names = FALSE)
write.csv(posterior.samples$muK,paste(runDir,"posterior_samples_muK.csv", sep = ""), row.names = FALSE)

