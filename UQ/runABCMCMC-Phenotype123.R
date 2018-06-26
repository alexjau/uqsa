# Uncertainty Quantification: Running ABC-MCMC with copulas
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

source('ScoringFunction.R')
source('loadTargets.R')
source('runModel.R')
source('SSsolutions.R')
source('copulaFunctions.R')
source('ABCMCMCFunctions.R')
source('PreCalibration.R')
library(ks)
library(VineCopula)
library(MASS)
library(R.utils)
library(R.matlab)

## use if on a cluster, changes to be made in 
## copulaFunctions.R, preCalibration.R and in this file
#library(doMC) 
#registerDoMC(20) # example with 20 cores throughout
#getDoParWorkers()

# default parameters and values from model
parVal <- c(0.0021524, 0.0088968, 0.028311, 0.041046, 0.046, 0.046, 0.046, 
            0.046, 0.046, 0.0016667, 0.016667, 0.008125, 0.092857, 15836, 
            22959, 1124, 4646, 0.028, 70, 800, 60, 600, 0.021, 0.04, 0.02, 
            0.002, 0.0044, 0.021, 0.021, 0.021, 0.021, 52, 2000, 4000, 5000, 
            2272.7, 0.04, 0.011, 0.011, 0.011, 0.011, 0.02, 0.002, 0.011, 
            0.0044, 2000, 4000, 5000, 2272.7, 0.73, 50, 0.05, 10, 10)

parNames <- c("kf*CaM*Ca", "kf*CaM_Ca1*Ca", "kf*CaM_Ca2*Ca", "kf*CaM_Ca3*Ca", 
              "kf*CaM*PP2B", "kf*CaM_Ca1*PP2B", "kf*CaM_Ca2*PP2B", "kf*CaM_Ca3*PP2B", 
              "kf*CaM_Ca4*PP2B", "kf*PP2B_CaM*Ca", "kf*PP2B_CaM_Ca1*Ca", "kf*PP2B_CaM_Ca2*Ca", 
              "kf*PP2B_CaM_Ca3*Ca", "KD*CaM_Ca3*Ca", "KD*CaM_Ca2*Ca", "KD*CaM_Ca1*Ca", 
              "KD*CaM*Ca", "KD*CaM_Ca4*PP2B", "KD*PP2B_CaM_Ca3*Ca", "KD*PP2B_CaM_Ca2*Ca", 
              "KD*PP2B_CaM_Ca1*Ca", "KD*PP2B_CaM*Ca", "kf*CaM*CaMKII", "kf*CaMKII_CaM_Ca3*Ca", 
              "kf*CaMKII_CaM_Ca2*Ca", "kf*CaMKII_CaM_Ca1*Ca", "kf*CaMKII_CaM*Ca", 
              "kf*CaM_Ca1*CaMKII", "kf*CaM_Ca2*CaMKII", "kf*CaM_Ca3*CaMKII", 
              "kf*CaM_Ca4*CaMKII", "KD*CaM_Ca4*CaMKII", "KD*CaMKII_CaM_Ca3*Ca", 
              "KD*CaMKII_CaM_Ca2*Ca", "KD*CaMKII_CaM_Ca1*Ca", "KD*CaMKII_CaM*Ca", 
              "kf*pCaMKII_Ca3*Ca", "kf*CaM*pCaMKIIaut", "kf*CaM_Ca1*pCaMKIIaut", 
              "kf*CaM_Ca2*pCaMKIIaut", "kf*CaM_Ca3*pCaMKIIaut", "kf*pCaMKII_Ca2*Ca", 
              "kf*pCaMKII_Ca1*Ca", "kf*CaM_Ca4*pCaMKIIaut", "kf*pCaMKII_Ca0*Ca", 
              "KD*pCaMKII_Ca3*Ca", "KD*pCaMKII_Ca2*Ca", "KD*pCaMKII_Ca1*Ca", 
              "KD*pCaMKII_Ca0*Ca", "KD*CaM_Ca4*pCaMKIIaut", "kautMax", "kf*PP1*pCaMKIIaut", 
              "kr*PP1*pCaMKIIaut", "kcat*PP1*pCaMKIIaut")

# parameters
parIdx <- c(14, 15, 16, 17, 18, 19, 20, 21, 22) 
parNames <- parNames[parIdx]

# scale to determine prior values
scale <- 1000

# prior
ll = parVal[parIdx]/scale
ul = parVal[parIdx]*scale
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale

# input values
np <- 26
xA <- seq(1,6,length.out = np) 
xB <- seq(0,7,length.out = np) 
xC <- seq(2,4,length.out = np)

# load targets
out <- loadTargets()
xtarget <- out$xtarget
ytarget <- out$ytarget
rm(out)

# no of samples
ns <- 1000 # no of samples required from each ABC-MCMC chain
npc <- 50000 # pre-calibration

# settings
p <- 0.01
nChains <- 20
delta <- 0.1

##############################
## Run MCMC for phenotype 1 ##
##############################
set.seed(7619201)
cat("Starting run for phenotype 1.\n")

# prior
ll = parVal[parIdx]/scale
ul = parVal[parIdx]*scale
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale

# start with independence copula for phenotype 1
out <- makeIndepCopula(ll, ul)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(10^xA, parIdx, xtarget[[1]], ytarget[[1]], npc, U, Z, copula, 'A')
sfactor <- 0.1 # scaling factor
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

drawsA <- lapply(1:nChains, function(i) ABCMCMC(10^xA, parIdx, ns, xtarget[[1]], ytarget[[1]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'A', ll, ul))

## on a cluster
#drawsA <- mclapply(1:nChains, function(i) ABCMCMC(10^xA, parIdx, ns, xtarget[[1]], ytarget[[1]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'A', ll, ul),mc.preschedule = FALSE, mc.cores = 20)

# put together
drawsA <- do.call("rbind", drawsA)
pick <- !apply(drawsA, 1, function(rw) all(rw==0))
drawsA <- drawsA[pick,]

## to save output for each step 
#draws <- drawsA
#outFile <- sprintf('Draws-Phenotype1-ABCMCMC-Scale%d.RData', scale)
#save(draws, parNames, file=outFile)

##############################
## Run MCMC for phenotype 2 ##
##############################
cat("Starting run for phenotype 2.\n")

# fit copula for draws already found for phenotype 1
out <- fitCopula(drawsA, ll, ul)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(10^xB, parIdx, xtarget[[2]], ytarget[[2]], npc, U, Z, copula, 'B')
sfactor <- 0.05 #scaling factor
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

drawsB <- lapply(1:nChains, function(i) ABCMCMC(10^xB, parIdx, ns, xtarget[[2]], ytarget[[2]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'B', ll, ul))

## on a cluster
#drawsB <- mclapply(1:nChains, function(i) ABCMCMC(10^xB, parIdx, ns, xtarget[[2]], ytarget[[2]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'B', ll, ul),mc.preschedule = FALSE, mc.cores = 20)

# put together
drawsB <- do.call("rbind", drawsB)
pick <- !apply(drawsB, 1, function(rw) all(rw==0))
drawsB <- drawsB[pick,]

# filter
pick <- apply(drawsB,1, function(u){out <- runModel(u, parIdx = parIdx, input = 10^xA, rInd = 'A');getMaxScore(xtarget[[1]], ytarget[[1]], out$xx, out$yy)}) <= delta
drawsB <- drawsB[pick,]  

## to save output for each step 
#draws <- drawsB
#outFile <- sprintf('Draws-Phenotype12-ABCMCMC-Scale%d.RData', scale)
#save(draws, parNames, file=outFile)

##############################
## Run MCMC for phenotype 3 ##
##############################
cat("Starting run for phenotype 3.\n")

# fit copula for draws already found for phenotype 1,2
out <- fitCopula(drawsB, ll, ul)
copula <- out$copula
U <- out$U
Z <- out$Z
Y <- out$Y

# precalibration all at once
out1 <- preCalibration(10^xC, parIdx, xtarget[[3]], ytarget[[3]], npc, U, Z, copula, 'C1')
sfactor <- 0.05 #scaling factor
out2 <- getMCMCPar(out1$prePar, out1$preDelta, p, sfactor, delta, nChains)
Sigma <- out2$Sigma
startPar <- out2$startPar

drawsC1 <- lapply(1:nChains, function(i) ABCMCMC(10^xC, parIdx, ns, xtarget[[3]], ytarget[[3]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'C1', ll, ul))

## on a cluster
#drawsC1 <- mclapply(1:nChains, function(i) ABCMCMC(10^xC, parIdx, ns, xtarget[[3]], ytarget[[3]], startPar[i,], Sigma, delta, U, Z, Y, copula, 'C1', ll, ul),mc.preschedule = FALSE, mc.cores = 20)

# put together
drawsC1 <- do.call("rbind", drawsC1)
pick <- !apply(drawsC1, 1, function(rw) all(rw==0))
drawsC1 <- drawsC1[pick,]

# filter
pick1 <- apply(drawsC1,1, function(u){out <- runModel(u, parIdx = parIdx, input = 10^xA, rInd = 'A');getMaxScore(xtarget[[1]], ytarget[[1]], out$xx, out$yy)}) <= delta

pick2 <- apply(drawsC1,1, function(u){out <- runModel(u, parIdx = parIdx, input = 10^xB, rInd = 'B');getMaxScore(xtarget[[2]], ytarget[[2]], out$xx, out$yy)}) <= delta

drawsC1 <- drawsC1[pick1 & pick2,]

draws <- drawsC1
outFile <- sprintf('Draws-Phenotype123-Scale%d.RData', scale)
save(draws, parNames, file=outFile)

outFile <- sprintf('Draws-Phenotype123-Scale%d.mat', scale)
writeMat(outFile, draws=draws, parNames=parNames)
