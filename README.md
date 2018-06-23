# uqsa
Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

The code distributed in this repository implements the methodology presented in the paper "Uncertainty quantification, propagation and
characterization by Bayesian analysis combined with global sensitivity analysis applied to dynamical intracellular pathway models" by Eriksson & Jauhiainen et al (2018). The code is distributed under the GNU General Public License v3.0.

The UQ folder contains R scripts to run the uncertainty quantification method (ABC-MCMC with copulas). The packages ks, VineCopula, MASS, and R.utils are required. The main script to run is called runABCMCMC-Phenotype123.R. 

