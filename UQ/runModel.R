# Uncertainty Quantification: Model simulation
# Copyright (C) 2018 
# Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# Olivia Eriksson
# Andrei Kramer
# Anu G Nair

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

runModel <- function(tpar, parIdx, input, rInd){
  
  parDefVal <- c(0.0021524, 0.0088968, 0.028311, 0.041046, 0.046, 0.046, 0.046, 
              0.046, 0.046, 0.0016667, 0.016667, 0.008125, 0.092857, 15836, 
              22959, 1124, 4646, 0.028, 70, 800, 60, 600, 0.021, 0.04, 0.02, 
              0.002, 0.0044, 0.021, 0.021, 0.021, 0.021, 52, 2000, 4000, 5000, 
              2272.7, 0.04, 0.011, 0.011, 0.011, 0.011, 0.02, 0.002, 0.011, 
              0.0044, 2000, 4000, 5000, 2272.7, 0.73, 50, 0.05, 10, 10)
  
  pars <- parDefVal
  pars[parIdx] <- 10^tpar
  yy <- numeric(length(input))
  
  if(rInd=='A'){
    yy <- SSsolutionFigA(pars, input) #param, Ca
  }else if(rInd=='B'){
    yy <- SSsolutionFigB(pars, input, 100) #param, CaM, PP2B
  }else if (rInd=='C1'){
    yy<- SSsolutionFigC(pars, input, 30,3) #param, Ca, totalCaM, totalPP2B
  }
  xx = log10(input)
  return(list(xx=xx, yy=yy))
}
