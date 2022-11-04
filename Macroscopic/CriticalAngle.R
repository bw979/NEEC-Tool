library(MASS)
library(car)
library(VGAM)
library(tidyverse)

### CONSTANTS ###
ab <- 5.2917721067E-11   #Bohr radius [m]
elec <- 1.60217662E-19;

### INPUTS ### 
A1 <- 84; Z1 <- 37; Z2 <- 14; EbeamMeVu <- 6; d <- 3.3
EbeamMev <- EbeamMeVu * 84 ; Ebeam <- EbeamMev * elec * 1E6
atf <- 0.885*ab*(Z1^0.5 + Z2^0.5)^(-2/3) #Thomas fermi screening radius [m]


### FUNCTION returns c(axial critical angle, planar critical angle)
# input units [beam A][beam Z1][target Z2][MeV/u][Angstrom]
crit_ang <- function(A1,Z1,Z2,Ebeam,d){
  ### Constants ###
  e <- 1.60217662E-19
  dmeters <- d*1E-10 
  Ejoules  <- Ebeam*A1*1E6*e
  
  ## Axial critical angle
  theta_c_axial <- ((2*Z1*Z2*e^2) / (Ejoules * dmeters))^0.5
  
  ## PLANAR
  atf <- 0.885*ab*(Z1^0.5 + Z2^0.5)^(-2/3) #Thomas fermi screening radius [m]
  Econd <- (2*Z1*Z2*d*elec^2)/(atf^2)
  npl <- 9.6E18  #atoms/m^2
  theta_c_pl <- ((2*pi*Z1*Z2*elec^2*npl*atf)/Ejoules)^0.5
  
  ## Output in rads
  return(c(theta_c_axial, theta_c_pl))
}


