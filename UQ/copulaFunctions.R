# Uncertainty Quantification: Copula functions for ABC-MCMC
# Copyright (C) 2018 Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

fitCopula <- function(X,ll,ul){
  
  ncx <- ncol(X)
  ns <- nrow(X)
  eps <- 0.1
  npoints <- 1000 
  #npoints <- 5000 # possible on a cluster
  
  # randomly pick sample points
  if(ns > npoints){
    I <- sample(1:ns, npoints, replace=FALSE)
  }else{
    I <- 1:ns
  }
  # add max and min
  I <- c(I, apply(X, 2, which.max)) 
  I <- c(I, apply(X, 2, which.min))
  I <- unique(I)
  
  Z <- U <- Y <-  matrix(NA, length(I), ncx)
  
  # must evaluate in real datapoints to 
  # keep connection between params
  # this is a normal kernel, looks similar 
  # to using the ecdf function
  for(i in 1:ncx){
    minx <- min(X[,i])
    maxx <- max(X[,i])
    ls <- max(ll[i],minx-eps) 
    us <- min(ul[i],maxx+eps)
    U[,i] <- X[I,i]
    Z[,i] = kcde(X[,i], xmin=ls, xmax=us, eval.points = X[I,i])$estimate
    Y[,i] = kde(X[,i], xmin=ls, xmax=us, eval.points = X[I,i])$estimate
  }
  
  # fit copula
  vineCop <- RVineStructureSelect(Z,indeptest = T)
  # vineCop <- RVineStructureSelect(Z,indeptest = T, mc.cores=20) # on a cluster
  return(list(copula=vineCop, U=U, Z=Z, Y=Y))
}


makeIndepCopula <- function(ll, ul){
  npoints <- 1000
  np <- length(ll)
  Z <- U <- Y <- matrix(NA, npoints, np)
   for(i in 1:np){
    minx <- ll[i]
    maxx <- ul[i]
    U[,i] <- seq(minx, maxx, length.out = npoints)
    Z[,i] <- seq(0,1, length.out=npoints) 
    Y[,i] <- rep(1/(maxx-minx), npoints) 
  }
  vineCop <- RVineStructureSelect(Z, family=0)
  return(list(copula=vineCop, U=U, Z=Z, Y=Y))
}

