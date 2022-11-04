# source("PDB/PDB_ResonanceStrength_Calculator.R")
# source("PDB/PDB_WeisskopfCalculator.R")

###CASE###
# #case_file<-"Theory_PDB_GS.csv"
# #57Fe using IC_K/2
# Jd <- 1.5; Jat <- 0.5; Ji <- 0.5; IC <- 7.35/2; t0 <- (98.3E-9)/log(2)
# Etrans_keV <- 14.4129; Ebind_keV <- 5.135 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# A <- 57
# Ar_1 <- compute_Ar("M1", A, Etrans_keV, 0.0078)
# Ar_2 <- 0
# Ar <- Ar_1 + Ar_2

# #84Rb using IC_M/10
# Jd <- 6; Jat <- 1.5; Ji <- 5; IC <- 50.5/10; t0 <- 9E-9
# Etrans_keV <- 3.050; Ebind_keV <- 0.5 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# Ar <- compute_Ar("M1", 84, Etrans_keV, 0.08)
# 
# #173Yb using IC_K/2
# Jd <- 7/2; Jat <- 0.5; Ji <- 5/2; IC <- 5.63/2; t0 <- 46E-12
# Etrans_keV <- 78.647; Ebind_keV <- 71.5748 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# Ar <- compute_Ar("M1", 173, Etrans_keV, 0.117)
# compute_S(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV)
#T_PDB <- read_csv(case_file)
#T_PDB$T0_d <- as.numeric(as.character(T_PDB$T0_d))

### Run over input list and allocate values
# for(i in 1:length(T_PDB$Spdb_P)) {
#   if(T_PDB$compute[i] == TRUE) {
#   Sexp <- compute_S(T_PDB$Jd[i], T_PDB$Jat[i], T_PDB$Ji[i], T_PDB$IC_P[i], T_PDB$Ar_SI[i], T_PDB$Etrans[i], T_PDB$Ebind[i]) * 1E-24
#   #Gamma_neec <- hbar/(T_PDB$T0_d[i]*e)
#   Gamma_neec <- Gamma[i]


##Running through a T_PDB table ... Form table first with PDB_Trends.R
#T_PDB <-


######### EXPERIMENTAL INPUTS ########
# ### COLLISION DENSITY... using typical values from NIST EBIT:
# ## Using [palffy thesis] upper estimates  
# # interaction volume
 # V <- 1#cm^3
 # nion <- 3E5 #cm^-3
 # nelec <- 1E13 #cm^-3



# ##### SPECIES INPUT #####
# Jd <- 1.5; Jat <- 0.5; Ji <- 0.5; IC <- 7.35/2; t0 <- (98.3E-9)/log(2)
# Etrans_keV <- 14.4129; Ebind_keV <- 5.135 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# A <- 57
# Ar_1 <- compute_Ar("M1", A, Etrans_keV, 0.0078)
# Ar_2 <- 0
# Ar <- Ar_1 + Ar_2


# #### upper limit F(E)
# s <- 50
# Ee <- 100000000
# Gamma_neec_eV <- 1E-5
# f_gam <- dnorm(Ee, mean=Ee, sd=s)
# #probability of a resonant collision
# f_gam * Gamma_neec_eV

# ## s is spread in electron beam energy
# length_trap <- 2 #cm
# R_elec <- 37E-4 #cm 
# V <- (R_elec^2)*pi*length_trap
# nion <- 1E9 #cm^-3
# nelec <- 1E11 #cm^-3
# s <- 50 #eV


                            #beV keV cm^3 cm^-3 cm^-3 eV
compute_EBIT_RATE <- function(S, Ee, Gamma_neec_eV, OUTPUT) {

  ## s is spread in electron beam energy
  length_trap <- 2 #cm
  R_elec <- 37E-4 #cm 
  V <- (R_elec^2)*pi*length_trap
  nion <- 1E9 #cm^-3
  nelec <- 1E11 #cm^-3
  s <- 50 #eV  
  
ncoll <- nion * nelec * V


### Dont need the width at the moment ###
#t0 <- (Thalf)/log(2)
#Etrans_keV <- 14.4129; Ebind_keV <- 5.135 #from nist
#Gamma_neec <- hbar/(t0*e)
#A <- 57
# Ar_1 <- compute_Ar("M1", A, Etrans_keV, 0.0078)
# Ar_2 <- 0
# Ar <- Ar_1 + Ar_2
 




### RESONANCE VELOCITY ###
vres_rel <- sqrt(1 - (1/((1+(Ee*1E3) / me))^2))  #c
vres <- 100*c*vres_rel #cm s^-1
# vres_rel <- sqrt(1 - (1/((1+(Wins_all$Ee*1E3) / me))^2))  #c
# vres <- 100*c*vres_rel #cm s^-1
#compute the resonance strength
#compute_S <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){
##
#Sexp_beV <- compute_S(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV) #beV
Sexp <- S * 1E-24  #cm^2 eV
##Lorentzian
Lor_max <- 2/(pi * Gamma_neec_eV)

### DETECTABLE FRACTION ###
#fraction that go into detectors
fp <- 1
#detector efficiency
ef <-  1 # 2E-3

## Gaussian Energy electron beam ##
#s<- 50  #...spread in electron beam energy in eV
#value of the PDF at the resonance
f_gam <- dnorm(Ee, mean=Ee, sd=s)
#probability of a resonant collision
Probability <- f_gam * Gamma_neec_eV
#proper probability
#proper_Probability <- pnorm(Ee + (Gamma_neec / 2), mean = Ee, sd = s) - pnorm(Ee - (Gamma_neec / 2), mean = Ee, sd = s)

#### How does the electron energy PDF combine with the ion beam energy PDF

rneec_gauss <-  ncoll * vres * Sexp * Lor_max * Probability *fp*ef
#using the proper integration
#Rneec_compare <- ncoll*vres*Sexp* Lor_max * proper_Probability *fp *ef
# using the algebraic simplification
#rneec_gauss <- ncoll * vres * Sexp * (sqrt(2)* pi^(3/2))/ s 
#rneec_gauss_Stot <-

Nneec_pw <- rneec_gauss * 60*60*24*7


### Will only NEED p as one only aims at a single resonance channel
if(OUTPUT == "Rate"){
  return(rneec_gauss)
}
# if(OUTPUT == "tot"){
#   return(rneec_gauss_Stot)
# }
if(OUTPUT == "N_neec_pw"){
  return(Nneec_pw) 
}
}




## Part of the T_PDB table for loop 
#}
#}







