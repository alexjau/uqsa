# Uncertainty Quantification: Model reductions
# Copyright (C) 2018 
# Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
# Olivia Eriksson
# Sara Maad Sasane
# Carolina Sartorius

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


SSsolutionFigA <- function(pars,Ca){

  # Free parameters
  KD1 = pars[17]; #[KD*CaM*Ca] = KD1
  KD2 = pars[16]; #[KD*CaM_Ca1*Ca] = KD2
  KD3 = pars[15]; #[KD*CaM_Ca2*Ca] = KD3
  KD4 = pars[14]; #[KD*CaM_Ca3*Ca] = KD4
  
  MolCaPerMolCaM <- (Ca*(4*Ca^3 + 3*KD4*Ca^2 + 2*KD3*KD4*Ca+ KD2*KD3*KD4))/(Ca^4 + KD4*Ca^3 + KD3*KD4*Ca^2 + KD2*KD3*KD4*Ca + KD1*KD2*KD3*KD4);
 
  return(MolCaPerMolCaM)
}

SSsolutionFigB <- function(pars,totalCaM,totalPP2B){
  
  # Free parameters
  KD1 = pars[17]; #[KD*CaM*Ca] = KD1
  KD2 = pars[16]; #[KD*CaM_Ca1*Ca] = KD2
  KD3 = pars[15]; #[KD*CaM_Ca2*Ca] = KD3
  KD4 = pars[14]; #[KD*CaM_Ca3*Ca] = KD4
  KD9 = pars[18]; #[KD*CaM_Ca4*PP2B] = KD9
  
  # Parameters defined by rules
  KD8 = (KD4 * KD9)/pars[19]; #[KD*CaM_Ca3*PP2B] = KD8
  KD7 = (KD3 * KD8)/pars[20]; #[KD*CaM_Ca2*PP2B] = KD7
  KD6 = (KD2 * KD7)/pars[21]; #[KD*CaM_Ca1*PP2B] = KD6
  KD5 = (KD1 * KD6)/pars[22]; #[KD*CaM*PP2B] = KD5
  
  MolCaMPerMolPP2B <- -(totalCaM*(KD5 - totalPP2B + totalCaM - (KD5^2 + 2*KD5*totalPP2B + 2*KD5*totalCaM + totalPP2B^2- 2*totalPP2B*totalCaM + totalCaM^2)^(1/2)))/(totalPP2B*(KD5 + totalPP2B - totalCaM + (KD5^2 + 2*KD5*totalPP2B + 2*KD5*totalCaM + totalPP2B^2 - 2*totalPP2B*totalCaM + totalCaM^2)^(1/2)))
  
  return(MolCaMPerMolPP2B)
}


SSsolutionFigC<-function(pars, Ca, totalCaM, totalPP2B){
  
  
  # Free parameters
  KD1 = pars[17]; #[KD*CaM*Ca] = KD1
  KD2 = pars[16]; #[KD*CaM_Ca1*Ca] = KD2
  KD3 = pars[15]; #[KD*CaM_Ca2*Ca] = KD3
  KD4 = pars[14]; #[KD*CaM_Ca3*Ca] = KD4
  KD9 = pars[18]; #[KD*CaM_Ca4*PP2B] = KD9
  
  
  # Parameters defined by rules
  KD8 = (KD4 * KD9)/pars[19]; #[KD*CaM_Ca3*PP2B] = KD8
  KD7 = (KD3 * KD8)/pars[20]; #[KD*CaM_Ca2*PP2B] = KD7
  KD6 = (KD2 * KD7)/pars[21]; #[KD*CaM_Ca1*PP2B] = KD6
  KD5 = (KD1 * KD6)/pars[22]; #[KD*CaM*PP2B] = KD5
  

    m1_7 <- (2*Ca*(Ca^4*KD5*KD6*KD7*KD8 + Ca^3*KD4*KD5*KD6*KD7*KD9  - Ca^4*KD5*KD6*KD7*KD8*KD9 + Ca^2*KD3*KD4*KD5*KD6*KD8*KD9  - Ca^3*KD4*KD5*KD6*KD7*KD8*KD9 - Ca^2*KD3*KD4*KD5*KD6*KD7*KD8*KD9 + Ca*KD2*KD3*KD4*KD5*KD7*KD8*KD9 + KD1*KD2*KD3*KD4*KD6*KD7*KD8*KD9 - KD1*KD2*KD3*KD4*KD5*KD6*KD7*KD8*KD9 - Ca*KD2*KD3*KD4*KD5*KD6*KD7*KD8*KD9)) /(Ca^2*KD5*KD6*KD7*KD8 + Ca^3*KD5*KD6*KD7*KD8 - Ca^4*KD5*KD6*KD7*KD8 - 2*Ca^3*KD4*KD5*KD6*KD7*KD9 + Ca^4*KD5*KD6*KD7*KD8*KD9 - 2*Ca^2*KD3*KD4*KD5*KD6*KD8*KD9 + Ca^3*KD4*KD5*KD6*KD7*KD8*KD9 + Ca^2*KD3*KD4*KD5*KD6*KD7*KD8*KD9 - 2*Ca*KD2*KD3*KD4*KD5*KD7*KD8*KD9 - 2*KD1*KD2*KD3*KD4*KD6*KD7*KD8*KD9 + KD1*KD2*KD3*KD4*KD5*KD6*KD7*KD8*KD9 + Ca*KD2*KD3*KD4*KD5*KD6*KD7*KD8*KD9)
    
    m2_7 <-  -(Ca^2*KD5*KD6*KD7*KD8 + Ca^3*KD5*KD6*KD7*KD8 - Ca^4*KD5*KD6*KD7*KD8 - 2*Ca^3*KD4*KD5*KD6*KD7*KD9 + Ca^4*KD5*KD6*KD7*KD8*KD9 - 2*Ca^2*KD3*KD4*KD5*KD6*KD8*KD9 + Ca^3*KD4*KD5*KD6*KD7*KD8*KD9 + Ca^2*KD3*KD4*KD5*KD6*KD7*KD8*KD9 - 2*Ca*KD2*KD3*KD4*KD5*KD7*KD8*KD9 - 2*KD1*KD2*KD3*KD4*KD6*KD7*KD8*KD9 + KD1*KD2*KD3*KD4*KD5*KD6*KD7*KD8*KD9 + Ca*KD2*KD3*KD4*KD5*KD6*KD7*KD8*KD9)/(Ca*KD5*KD6*KD7*KD8)
    
    ActivePP2BPercentage <- (200*(Ca^4/(2*Ca^4 + 2*Ca^3 + 2*Ca^2 + 2*m2_7*Ca + m1_7*m2_7) - (((2*Ca^4)/(2*Ca^4 + 2*Ca^3 + 2*Ca^2 + 2*m2_7*Ca + m1_7*m2_7) + (Ca^4*(totalPP2B - totalCaM))/(Ca^4 + Ca^3 + Ca^2 + m2_7*Ca + m1_7*m2_7))^2 + (8*Ca^8*totalCaM)/((Ca^4 + Ca^3 + Ca^2 + m2_7*Ca + m1_7*m2_7) *(2*Ca^4 + 2*Ca^3 + 2*Ca^2 + 2*m2_7*Ca + m1_7*m2_7)))^(1/2)/2 + (Ca^4*(totalPP2B - totalCaM))/(2*(Ca^4 + Ca^3 + Ca^2 + m2_7*Ca + m1_7*m2_7)))) /((((2*Ca^4)/(2*Ca^4 + 2*Ca^3 + 2*Ca^2 + 2*m2_7*Ca + m1_7*m2_7) - (((2*Ca^4)/(2*Ca^4 + 2*Ca^3 + 2*Ca^2 + 2*m2_7*Ca + m1_7*m2_7) + (Ca^4*(totalPP2B - totalCaM))/(Ca^4 + Ca^3 + Ca^2 + m2_7*Ca + m1_7*m2_7))^2 + (8*Ca^8*totalCaM)/((Ca^4 + Ca^3 + Ca^2 + m2_7*Ca + m1_7*m2_7) *(2*Ca^4 + 2*Ca^3 + 2*Ca^2 + 2*m2_7*Ca + m1_7*m2_7)))^(1/2) + (Ca^4*(totalPP2B - totalCaM))/(Ca^4 + Ca^3 + Ca^2 + m2_7*Ca + m1_7*m2_7)) *(2*Ca^4 + 2*Ca^3 + 2*Ca^2 + 2*m2_7*Ca + m1_7*m2_7))/(2*Ca^4) - 2)
  
  # KD-parameters 
  #[KD*CaM_Ca3*Ca] = KD4
  #[KD*CaM_Ca2*Ca] = KD3
  #[KD*CaM_Ca1*Ca] = KD2
  #[KD*CaM*Ca] = KD1
  #[KD*CaM_Ca4*PP2B] = KD9
  #[KD*CaM_Ca3*PP2B] = KD8
  #[KD*CaM_Ca2*PP2B] = KD7
  #[KD*CaM_Ca1*PP2B] = KD6
  #[KD*CaM*PP2B] = KD5
  return(ActivePP2BPercentage)
}

