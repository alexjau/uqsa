# Uncertainty Quantification: Precalibration for ABC-MCMC
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

preCalibration <- function(input, parIdx, xtarget, ytarget, npc, U, Z, copula, rInd){
  
  np <- length(parIdx)
  
  R <- RVineSim(npc, copula)
  prePar <- matrix(0, npc, np)
  for(i in 1:np){
    prePar[,i] = spline(Z[,i],U[,i],xout=R[,i])$y
  }
	
	fun <- function(i) {
		tpar <- prePar[i,]
    	invisible(capture.output(out <- withTimeout(runModel(tpar, parIdx, input, rInd), timeout=25, onTimeout="silent")))
		tmp <- ifelse(is.null(out), NA, getMaxScore(xtarget, ytarget, out$xx, out$yy))
		return(tmp)
	}

  preDelta <- lapply(1:npc, fun)
  # preDelta <- mclapply(1:npc, fun, mc.preschedule = FALSE, mc.cores = 20) # on a cluster
  preDelta <- unlist(preDelta)

  return(list(preDelta=preDelta, prePar=prePar))
}



getMCMCPar <- function(prePar, preDelta, p, sfactor, delta, nChains){
  prePar <- prePar[!is.na(preDelta),]
  preDelta <- preDelta[!is.na(preDelta)]
  nk <- nrow(prePar)*p
  pick1  <- grep(TRUE, (preDelta <= delta)) # pick all pars that meet threshold
  pick2 <- order(preDelta, decreasing = F)[1:nk] # pick top p percent
  
  if(length(pick1)>length(pick2)){
    pick <- pick1
  }else{
    pick <- pick2
  }

  Scorr <- cor(prePar[pick,])*0.8 #tone down corrs
  diag(Scorr) <- 1
  sdv <- apply(prePar[pick,], 2, sd)
  Sigma <- sfactor * Scorr * tcrossprod(sdv)
 
  startPar <- prePar[sample(pick, nChains, replace = FALSE),]
  list(Sigma=Sigma, startPar=startPar)
}
