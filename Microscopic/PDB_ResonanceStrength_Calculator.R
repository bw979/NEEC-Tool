library(tidyverse)
library(stringr)
library(plotly)
#for the dfactorial function
#library(phangorn)



####### Calculates S and computes the NEEC rate based on this ########
####### ATOMIC UNITS within the S function, SI units elsehere using eV where appropriate) ##############
##output resonance strength

## Fundamental Constants
#AU to SI
e <- 1.602E-19 #C
me <- 0.51099895000E6 #eV/c^2
hbar <- 1.054571817E-34	#J.s
Eh <- 4.359744722E-18 ### one hartree is that many Joules
Eh_eV <- 27.211386 ## one hartree is that many electron volts

#SI to a.u conversion
time_conv <- 2.418884326E-17 # this many seconds is an atomic unit ...hbar/Eh
length_conv <- 5.291772E-11 # this many m is one atomic unit ...length of one Bohr radius a0
barn <- 1E-28 #m^2
me.mp <- 1/(1836.15267343)
mp.me <- 1/me.mp
c <- 299792458 #ms^-1



##Ar - Radiative transition rate (first read B(lam) from GAMMA continuation card)
#RETURNS W-Eckart estimate function (SI units atm)
#INPUT: L, B(multiples of W.U.), Etrans (eV)
  # compute_Ar <- function(L, B, Etrans_eV){
  #   Etrans_MeV <- Etrans_eV*(10^-6)
  #   c <- 137.036
  #   Ar <- ((8*pi*(L+1))/(L*((dfactorial(2*L+1))^2)))*((Etrans_MeV^(2*L+1))/(c))*B
  #   return(Ar)
  # }

#computes S with input Everything in a.u. except Etrans_keV and Ar in s^-1
compute_S <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){
  if(Jat == 0) {Jat = 0.5}
  # obs$Jd, obs$Ji, obs$Jat, , obs$Ar
  ## The transition energy into AU
  Etrans_J <- 1E3*Etrans_keV * e
  Etrans_AU <- Etrans_J/Eh
  ## The continuum electron energy into AU
  Ebind_J <- 1E3*Ebind_keV *e
  Ebind_AU <- Ebind_J/Eh
  Econt_AU <- Etrans_AU - Ebind_AU
  # return(Etrans_AU)
  c <- 137.036
  p <- sqrt( ((Econt_AU + c^2)^2) - (c^4) ) / c
  Yneec_SI <-  ( (2*Jd + 1)*(2* Jat + 1) / (2*(2*Ji + 1)) ) * IC * Ar
  Yneec_AU <- Yneec_SI * time_conv
  Spdb_AU <- ( (2*(pi^2) ) / (p^2) ) * Yneec_AU
  Spdb_beV <- Spdb_AU*((length_conv^2)*Eh_eV)*(1/barn)
  return(Spdb_beV)
}

compute_S_scaled <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV, Vi_q0_keV, n, Element){
  if(Jat == 0) {Jat = 0.5}
  # obs$Jd, obs$Ji, obs$Jat, , obs$Ar
  ## The transition energy into AU
  Etrans_J <- 1E3*Etrans_keV * e
  Etrans_AU <- Etrans_J/Eh
  ## The continuum electron energy into AU
  Ebind_J <- 1E3*Ebind_keV *e
  Ebind_AU <- Ebind_J/Eh
  Econt_AU <- Etrans_AU - Ebind_AU
  # return(Etrans_AU)
  c <- 137.036
  p <- sqrt( ((Econt_AU + c^2)^2) - (c^4) ) / c
  
  #read the subshell that should be used for that atom
  
  #scale the IC coeff
  IC_scaled <- (Ebind_keV/Vi_q0_keV)*(nhole/nmax)*IC
  
  Yneec_SI <-  ( (2*Jd + 1)*(2* Jat + 1) / (2*(2*Ji + 1)) ) * IC_scaled * Ar
  Yneec_AU <- Yneec_SI * time_conv
  Spdb_AU <- ( (2*(pi^2) ) / (p^2) ) * Yneec_AU
  Spdb_beV <- Spdb_AU*((length_conv^2)*Eh_eV)*(1/barn)
  return(Spdb_beV)
}

#computes S with input Everything in a.u. except Etrans_keV and Ar in s^-1
compute_S_check <- function(Jd, Jat,Ji,IC, Yn, Etrans_keV, Ebind_keV){
  if(Jat == 0) {Jat = 0.5}
  # obs$Jd, obs$Ji, obs$Jat, , obs$Ar
  ## The transition energy into AU
  Etrans_J <- 1E3*Etrans_keV * e
  Etrans_AU <- Etrans_J/Eh
  ## The continuum electron energy into AU
  Ebind_J <- 1E3*Ebind_keV *e
  Ebind_AU <- Ebind_J/Eh
  Econt_AU <- Etrans_AU - Ebind_AU
  # return(Etrans_AU)
  c <- 137
  p <- Econt_AU / c
  #Yneec_SI <-  ( (2*Jd + 1)*(2* Jat + 1) / (2*(2*Ji + 1)) ) * IC * Ar
  Yneec_AU <- Yn * time_conv
  Spdb_AU <- ( (2*(pi^2) ) / (p^2) ) * Yneec_AU
  Spdb_beV <- Spdb_AU*((length_conv^2)*Eh_eV)*(1/barn)
  return(Spdb_beV)
}

compute_S_check_Y <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){
  if(Jat == 0){Jat = 0.5}
  # obs$Jd, obs$Ji, obs$Jat, , obs$Ar
  ## The transition energy into AU
  Etrans_J <- 1E3*Etrans_keV * e
  Etrans_AU <- Etrans_J/Eh
  ## The continuum electron energy into AU
  Ebind_J <- 1E3*Ebind_keV *e
  Ebind_AU <- Ebind_J/Eh
  Econt_AU <- Etrans_AU - Ebind_AU
  # return(Etrans_AU)
  c <- 137.036
  p <- sqrt( ((Econt_AU + c^2)^2) - (c^4) ) / c
  Yneec_SI <-  ( (2*Jd + 1)*(2* Jat + 1) / (2*(2*Ji + 1)) ) * IC * Ar
  Yneec_AU <- Yneec_SI * time_conv
  Spdb_AU <- ( (2*(pi^2) ) / (p^2) ) * Yneec_AU
  Spdb_beV <- Spdb_AU*((length_conv^2)*Eh_eV)*(1/barn)
  return(Yneec_SI)
}


# compute_S_alphaTot <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){
#   if(Jat == 0) {Jat = 0.5}
#   # obs$Jd, obs$Ji, obs$Jat, , obs$Ar
#   ## The transition energy into AU
#   Etrans_J <- 1E3*Etrans_keV * e
#   Etrans_AU <- Etrans_J/Eh
#   ## The continuum electron energy into AU
#   Ebind_J <- 1E3*Ebind_keV *e
#   Ebind_AU <- Ebind_J/Eh
#   Econt_AU <- Etrans_AU - Ebind_AU
#   # return(Etrans_AU)
#   c <- 137.036
#   p <- sqrt( ((Econt_AU + c^2)^2) - (c^4) ) / c
#   Yneec_SI <-  ( (2*Jd + 1)*(2* Jat + 1) / (2*(2*Ji + 1)) ) * IC * Ar
#   Yneec_AU <- Yneec_SI * time_conv
#   Spdb_AU <- ( (2*(pi^2) ) / (p^2) ) * Yneec_AU
#   Spdb_beV <- Spdb_AU*((length_conv^2)*Eh_eV)*(1/barn)
#   return(Spdb_beV)
# }



## Mulitpolarity Selector(See Eexc search)
#find the right B to use in the Ar calc   

## Finish proper Ar rate calculator


## Fill in missing IC values with BRICC

# ##84Rb isomer case
# Eres_keV <- 3.050
# Ebind_keV <- 0.0565
# BW <-0.08
# Ji <-6
# Jd <-5
# Jat <-0.5
# IC <- 0.32
# Eres_MeV <- Eres_keV*(10^-3)
# #Ar_TM1 <- 1.779E13*
# B_SI <- 1.790 * BW
# Ar <- 1.779E13*(Eres_MeV^3)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# Sexp <- compute_S(Jd,Jat,Ji,IC,Ar,Eres_keV, Ebind_keV)
# ### OUTPUT : 0.01066692 ... Agrees with same calc done by Palffy




