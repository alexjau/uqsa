# uqsa
Uncertainty Quantification (UQ) and Sensitivity Analysis (SA)

The code distributed in this repository implements the methodology presented in the paper "Uncertainty quantification, propagation and
characterization by Bayesian analysis combined with global sensitivity analysis applied to dynamical intracellular pathway models" by Eriksson & Jauhiainen et al (2018). The code is distributed under the GNU General Public License v3.0.

The UQ folder contains R scripts to run the uncertainty quantification method (ABC-MCMC with copulas). The packages ks, VineCopula, MASS, R.utils, and R.matlab are required (the last package to save the output data also into a MATLAB compatible file). The main script to run is called runABCMCMC-Phenotype123.R. This script will fit the model to phenotypes 1-3 (as described in the paper), which we use as an illustrative test case. The resulting data is also uploaded in the folder, both in R and MATLAB format. We use phenotype 4 as our prediction dataset to illustrate the SA methodology. 

Th SA folder contains MATLAB scripts to run global sensitivity analysis. A version of MATLAB later than 2014a is required. The main script to run is called get_predictions_do_SA.m and requires access to the file Draws-Phenotype123-Scale1000.mat which is available in the UQ folder. 


