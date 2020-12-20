# pOXA48

This is the repository of scripts and data necessary to generate the figures in the main text and supplementary information of the article *"The distribution of plasmid fitness effects explains plasmid persistence in bacterial communities"*.  Experimental data for each strain in our collection can be found at *pOXA48/data/ODdata_strains/*.

The following scripts can be found in *pOXA48/code/m-files/*:

**Figure 5a-b: Modelling pOXA-48 fitness effects**
* Plotting script: plotFigure5ab.m

**Figure 5c: Fraction of plasmid-bearing as a function of conjugation rate**
* Simulation script: runFigure5c.m
* Plotting script: plotFigure5c.m

**Figure 5d: Fraction of plasmid-bearing cells as a function of the rate of HGT**
* Simulation script: runFigure5d.m
* Plotting script: plotFigure5d.m

**Figure 6: Modelling plasmid persistence in polymicrobial communities**
* Simulation script: runFigure6.m
* Plotting script: plotFigure6.m

**Figure 6e: Mean fraction of plasmid-bearing cells as a function of the number of strains in the community**
* Simulation script: runFigure6e.m
* Plotting script: plotFigure6e.m

**Supp Figure 7:In silico competition experiments**
* Plotting script: plotSuppFigure7.m

**Supplementary Figure 8.  Comparison between in silico and experimental competition experiments**
* Plotting script: plotSuppFigure8.m

**Supp Figure 9. Plasmid stability as a function of the conjugation rate and plasmid cost**
* Plotting script: plotSuppFigure9.m

**Supp Figure 10: Effect of conjugation and community complexity in plasmid population dynamics**
* Simulation script: runSuppFigure10.m
* Plotting script: plotSuppFigure10.m

**Supp Figure 15: Sampling the distribution of plasmid fitness effects**
* Simulation script: plotSuppFigure15.m

## Model parametrization

Parameters of the population dynamics model (specific affinity and cell efficiency) were jointly determined using a Markov chain Monte Carlo method with a Metropolis-Hastings sampler from growth curves measured as bacterial optical densities.  It also implements a data cloning algorithm to assess the identifiability of the parameters.

R scripts to perform the parametrization for the strain collection presented in this manuscript can be run by executing the script *runMe* that can be found in *pOXA48/code/MCMC/R/*

MCMC chains are exported into *.csv* files that can be found in *pOXA48/data/MCMC_chains/*.  

Diagnostic plots presented in Supplementary File 1 are produced with the following Matlab script:

**Supp File 1. MCMC Diagnostic plots**
* Plotting script: plotSuppFile1.m
* Gelman-Rubin statistics: runGelmanRubin.m
