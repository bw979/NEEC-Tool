# source("PDB/PDB_ResonanceStrength_Calculator.R")
# source("PDB/PDB_WeisskopfCalculator.R")

###CASE###
##Running through a T_PDB table ... Form table first with PDB_Trends.R
#T_PDB <-

# #57Fe using IC_K/2
# Jd <- 1.5; Jat <- 0.5; Ji <- 0.5; IC <- 7.35/2; t0 <- 98.3E-9
# Etrans_keV <- 14.4129; Ebind_keV <- 5.135 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# Ar <- compute_Ar("M1", 57, Etrans_keV, 0.0078)

# #57Fe using IC_K/2
# Jd <- 1.5; Jat <- 0.5; Ji <- 0.5; IC <- 7.35/2; t0 <- (98.3E-9)/log(2)
# Etrans_keV <- 14.4129; Ebind_keV <- 5.135 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# #A <- 57
#Ar_1 <- compute_Ar("M1", A, Etrans_keV, 0.0078)
#Ar_2 <- 0
#Ar <- Ar_1 + Ar_2





# #84Rb using IC_M/10
# Jd <- 6; Jat <- 1.5; Ji <- 5; IC <- 50.5/10; t0 <- 9E-9
# Etrans_keV <- 3.050; Ebind_keV <- 0.5 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# Ar <- compute_Ar("M1", 84, Etrans_keV, 0.08)

# #173Yb using IC_K/2
# Jd <- 7/2; Jat <- 0.5; Ji <- 5/2; IC <- 5.63/2; t0 <- 46E-12
# Etrans_keV <- 78.647; Ebind_keV <- 71.5748 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# Ar <- compute_Ar("M1", 173, Etrans_keV, 0.117)
# compute_S(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV)

# ## 129mSb using tabulated IC_L/8
# Jd <- 15/2; Jat <- 0.5; Ji <- 19/2; IC <- 2.72E4; t0 <- (2.2E-6)/log(2) # the halflife converted to mean lifetime
# Etrans_keV <- 1861.06 - 1851.31; Ebind_keV <- 8.670 #from nist
# Econt <- Etrans_keV - Ebind_keV;  Gamma_neec <- hbar/(t0*e)
# A <- 129; mion.me <- A*mp.me
# Ar_1 <- compute_Ar("E2", A, Etrans_keV, 1.96)
# Ar_2 <- 0
# Ar <- Ar_1 + Ar_2



#for(i in 1:length(T_PDB$Spdb_P)) {
#  if(T_PDB$compute[i] == TRUE) {
#    Sexp <- compute_S(T_PDB$Jd[i], T_PDB$Jat[i], T_PDB$Ji[i], T_PDB$IC_P[i], T_PDB$Ar_SI[i], T_PDB$Etrans[i], T_PDB$Ebind[i]) * 1E-24
#    Gamma_neec <- hbar/(T_PDB$T0_d[i]*e)


### COLLISION DENSITY... GSI paramters:
# interaction volume
V <- 1#cm^3
nion <- 1E5 #cm^-3
nelec <- 1E10 #cm^-3
ncoll <- nion * nelec * V 
#channelled fraction


### RESONANCE VELOCITY ###
vres <- sqrt(2*(Econt*1E3)/(me)) #c
vres <- 100*c*vres #cm s^-1
#compute the resonance strength
#compute_S <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){

### RESONANCE STRENGTH ###
Sexp_beV <- compute_S(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV)  #beV
Sexp <- Sexp_beV * 1E-24 #cm^2 eV
#Lorentzian
Lor_max <- 2/(pi * Gamma_neec)
sigma_max <- Sexp * Lor_max #b


### Detectable fraction###
#fraction that go into detector
fp <- 1
#detector efficiency
ef <- 1

## Gaussian electron energy ##
sb <- 18 #keV
s <- (1/mion.me)*sb * 1000   #...spread in electron beam energy in eV
#value of the PDF at the resonance
f_gam <- dnorm(Econt, mean=Econt, sd=s)
## probability of a resonant collision given spread in electron energy and beam energy
#proper probability
Probability <- f_gam * Gamma_neec 
proper_Probability <- pnorm(Econt + (Gamma_neec / 2), mean = Econt, sd = s) - pnorm(Econt - (Gamma_neec / 2), mean = Econt, sd = s)

#### How does the electron energy PDF combine with the ion beam energy PDF

rneec_gauss <-  ncoll * vres * Sexp * Lor_max * Probability *fp*ef
#using the proper integration
Rneec_compare <- 100E4 *nion*nelec * 100E-4 * sigma_max * (Gamma_neec/s) / (1E-3)
rneec_gauss
Rneec_compare

 
  








