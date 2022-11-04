library(tidyverse)
library(stringr)
library(plotly)
#for the dfactorial function
library(phangorn)

###### WHAT DOES IT DO: Produces AnN0_S
### Calculates the resonance strength for a certain spin orbit shell
# - Provides a lookup table of Internal Conversion coefficients (CC) of a certain gamma ray
# - Compares transitions between levels to a relevent gamma ray energy (for the isomeric case)

#TODO
# Ebind column
# Ji and Ei column (for isomeric)

#Columns names: Ei: Initial nuclear energy, Ed: Upper nuclear energy, 

#functions
source("tidy_functions.R")

#### Inputs 
# some "output/AnNi" csv file which contains nuclear level entries E_init, E_final, Jpi and atomic capture orbital, J_at
# IC coefficients for the relevant shell (atm will use total IC coeff)

#if(!file exists) 
#source("Continuation_Cleaner.R")
#else ND_Gamma

##### using total IC presently ######

##Read Input Data
ND_Gamma <- read_csv("Dependencies/ND_Gamma_Cont.csv") %>%
  distinct()
A0N0 <- read_csv("Outputs/AnN0/A0N0_beam_001_298.csv") %>%
  distinct()
ND <- read_csv("Dependencies/ND_TidyJpi.csv") %>%
  distinct()


## Prototype for transition resonance strength data.frame AnN0_S 
#tibble( AX = character(1), Ei = double(1), Ed = double(1), IC = double(1), Ji = double(1), Jd = double(1) )
A0N0_S <- A0N0 %>% rename(Ed = E, Jd = J, Thf = Thalf) %>% 
  select( -IAXtf, -Gtf ) %>%
  mutate( IC = double(1), Jat = double(1) , Ebind = double(1), S = double(1) )
rm(A0N0)


## create unique GS array 
##(sep file?) Get list of ground state nuclides and fill Ji and Ei column 
GS <- ND %>% filter(E == 0.0)
GS <- tidy_Jpi(GS) %>% 
  select(-Z, -M, -El, -DThalf, -Tentative_Jpi, -Multiple_Jpi, -DE) %>%
  rename(Ei = E, Ji = J, Thi = Thalf) %>%
  distinct() %>%
  filter(!is.na(Thi))

# count(GS, AX) %>%
#   filter(n>1)
# around 170 with n>1 ... not all GS has unique AX key yet

#JOIN#
## Adds the Ji and Ei column by doing a left outer join
A0N0_S_GS <- left_join(A0N0_S, GS, by = "AX") %>%
 mutate(Etrans = Ed - Ei) 

##### By how much should we round in order to match keys?? ####
d <- 2
A0N0_S_GS$Etrans <- round(A0N0_S_GS$Etrans, d)
ND_Gamma$Egam <- round(ND_Gamma$Egam, d)

## (sep file?) Lookup and merge atomic J_at info
#which spin-orbit shell are we in 
eJconv <- read_csv("Dependencies/elec_Jconv.csv")
ebinds <- read_csv("Dependencies/EBind4.csv")
##SHOULD BE A FOR LOOP FOR ALL AnN0 files
n <- 1
A0N0_S_GS$Jat <- eJconv$Jat[n]
# for(i in 1:length(A0N0_S_GS$Ebind)){
#   A0N0_S_GS$Ebind[i] <- ebind$==a     
# }

#JOIN#
##Assign IC column
#Do [filtering join] to get A0N0_S$IC column read from ND_Gamma
A0N0_S_GAM <- left_join(A0N0_S_GS, ND_Gamma, by = c("AX", "Etrans" = "Egam") ) %>%
  mutate(S_tot = double(1), type = character(1), Ar_SI = double(1)) %>%
  select(AX, M, El, Z, Ji, Thi, Ei, Jd, Ed,Thf, Jat, Ebind, IC, CC, Etrans, type, Continuation, Braw, Btype, B, Ar_SI, S_tot, -M, -El, -Z ) 
  #rename(Ei = E, Ji = J, Thi = Thalf) 


## Fill from correct ebinds column
ecol <- n+2

#seems to feed in numeric(0) atm
for(i in 1:length(A0N0_S_GAM$Ebind)){
  A0N0_S_GAM$Ebind[i] <- ebinds[which(ebinds$Symb==A0N0_S_GAM$El[i]), ecol]
}

# eb <- function(obs) {
#   obs$Ebind <- Ebind[which(Ebind$Symb==obs$El), ecol]
# }

  




####### S column of AnN0_S ########
####### ATOMIC UNITS ##############
##output resonance strength

# ##Constants
# e <- 1.6E-19
# Eh <- 4.359744722E-18 ### one hartree is that many Joules
# Eh_eV <- 27.211386 ## one hartree is that many electron volts
# ##... in a.u
# time_conv <- 2.418884326E-17 # this many seconds is an atomic unit
# length_conv <- 5.291772E-11 # this many Kg in one atomic unit
# barn <- 1E-28



##Ar - Radiative transition rate (first read B(lam) from GAMMA continuation card)
#RETURNS W-Eckart estimate function (SI units atm)
#INPUT: L, B(multiples of W.U.), Etrans (eV)
# compute_Ar <- function(L, B, Etrans_eV){
#   Etrans_MeV <- Etrans_eV*(10^-6)
#   c <- 137.036
#   Ar <- ((8*pi*(L+1))/(L*((dfactorial(2*L+1))^2)))*((Etrans_MeV^(2*L+1))/(c))*B
#   return(Ar)
# }


# #computes S with input Everything in a.u. except Etrans_eV
# compute_S <- function(Jd, Jat,Ji,IC, Ar, Etrans_eV){
#   # obs$Jd, obs$Ji, obs$Jat, , obs$Ar
#   Etrans_J <- Etrans_eV *e
#   Etrans_AU <- Etrans_J/Eh
#   # return(Etrans_AU)
#   c <- 137
#   p <- Etrans_AU / c
#   Yneec_SI <-  ( (2*Jd + 1)*(2* Jat + 1) / (2*(2*Ji + 1)) ) * IC * Ar
#   Yneec_AU <- Yneec_SI * time_conv
#   Spdb_AU <- ( (2*(pi^2) ) / (p^2) ) * Yneec_AU
#   Spdb_beV <- Spdb_AU*((length_conv^2)*Eh_eV)*(1/barn)
#   return(Spdb_beV)
# 
#   }

## Mulitpolarity Selector(See Eexc search)
#find the right B to use in the Ar calc

## Finish proper Ar rate calculator
#set mass
  



## Fill in missing IC values with BRICC

## Fill in the PDB tot and K columns




# ############################# CASES #################################
# ##84Rb isomer case
# Eres_eV <- 3050
# BW <-0.08
# Ji <-6
# Jd <-5
# Jat <-0.5
# IC <- 0.32
# Eres_MeV <- Eres_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- 1.790 * BW
# Ar <- 1.779E13*(Eres_MeV^3)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,Eres_eV)
# ### OUTPUT : 0.01066692 ... Agrees with same calc done by Palffy
# 
# 
# #### SOME PALFFY HYDROGENIC CASES
# #165Ho case
# En_eV <- 94700
# BW <- 0.275
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- 1.790 * BW
# Ar <- 1.779E13*(En_MeV^3)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# IC <- 3.0
# compute_S(7/2,0.5,9/2,IC,Ar,En_eV)
# ### OUTPUT : 10.09123 beV ... Palffy Value:0.884 beV
# 
# 
# # ## NO B or IC VALUE #### #173Yb case
# # Eres_eV <- 78640
# # BW <- 0.275
# # SP_WeisskopfM1_SI <- 1.79
# # Eres_MeV <- Eres_eV*(10^-6)
# # #Ar_TM1 <- 1.779E13*
# # B_SI <- 1.790 * BW
# # Ar <- 1.779E13*(Eres_MeV^3)*B_SI # PER SECOND ... works spot on
# # #compute_Ar(1,B_SI,3200) #  Should match
# # #Ar_SI <- Ar / time_conv
# # IC <- 3.0
# # compute_S(7/2,0.5,5/2,IC,Ar,Eres_eV)
# # ### OUTPUT : 
# 
# # ## NO B or IC VALUE #### #185Re case
# # Eres_eV <- 
# # BW <- 0.275
# # SP_WeisskopfM1_SI <- 1.79
# # Eres_MeV <- Eres_eV*(10^-6)
# # #Ar_TM1 <- 1.779E13*
# # B_SI <- 1.790 * BW
# # Ar <- 1.779E13*(Eres_MeV^3)*B_SI # PER SECOND ... works spot on
# # #compute_Ar(1,B_SI,3200) #  Should match
# # #Ar_SI <- Ar / time_conv
# # IC <- 3.0
# # compute_S(7/2,0.5,5/2,IC,Ar,Eres_eV)
# # ### OUTPUT : 
# 
# #187Re case
# En_eV <- 134000.243
#   BW <-0.260
#   Ji <-5/2
#   Jd <-7/2
#   Jat <-0.5
#   IC <- 2.2
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- 1.790 * BW
# Ar <- 1.779E13*(En_MeV^3)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# ### OUTPUT : 16 beV ... Palffy Value:1.16 beV
# 
# 
# 
# #####  55Mn M1
# En_eV <- 125.95E3
# BW <-0.0417
# Ji <-5/2
# Jd <-7/2
# Jat <-0.5
# IC <- 0.01691
# Eres_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- 1.790 * BW
# Ar <- 1.779E13*(En_MeV^3)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# ### L shell IC OUTPUT : 0.01903 beV ... Palffy Value: 9.22E-4 beV
# ### Total IC OUTPUT :   0.01903 beV
# 
# 
# 
# 
# ##### E2 #####
# 
# #####  129Sb 17min isomer E2
# En_eV <- 9.76E3
# BW <-1.96
# Ji <-19/2
# Jd <-15/2
# Jat <-1.5
# IC <- 5.59E3
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 129^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# ### L shell OUTPUT : 0.01903 beV ... Palffy Value: 1.6E-3 beV
# ### M shell OUTPUT : 0.00586 beV ... Palffy Value: 3E-5 beV
# 
# 
# ## 164 Dy
# En_eV <- 73.392E3
# BW <- 211
# Ji <-0
# Jd <-2
# Jat <-0.5
# IC <- 8.89
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 164^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# 
# 
# ## 170Er
# En_eV <- 78.591E3
# BW <- 208
# Ji <-0
# Jd <-2
# Jat <-0.5
# IC <- 7.47
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 164^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# 
# ## 174Yb
# En_eV <- 76.471E3
# BW <- 201
# Ji <-0
# Jd <-2
# Jat <-0.5
# IC <- 9.43
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 164^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# 
# 
# ## 174Yb
# En_eV <- 76.471E3
# BW <- 201
# Ji <-0
# Jd <-2
# Jat <-0.5
# IC <- 9.43
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 164^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# 
# ## 154Gd
# En_eV <- 123.071E3
# BW <- 157
# Ji <-0
# Jd <-2
# Jat <-0.5
# IC <- 2.297
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 164^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# 
# ## 156Gd
# En_eV <- 88.966E3
# BW <- 189
# Ji <-0
# Jd <-2
# Jat <-1.5
# IC <- 0.24
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 164^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# 
# ## 62Dy
# En_eV <- 80.066E3
# BW <- 204
# Ji <-0
# Jd <-2
# Jat <-0.5
# IC <- 6.14
# En_MeV <- En_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- (5.94E-2)* 164^(4/3) * BW
# Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # #Data to plot
# # df=ebinds
# # 
# # #### EDIT HERE TO MAKE IT LOOK NICE ####
# # p <- plot_ly(df, x =~Z, y =~, z=~output, type="scatter3d", mode="markers+lines", line=list(width=2, dash="dash"), opacity=1.0, size=4, name=~El)
# # htmlwidgets::saveWidget(p, "EvenEven_Plots.html")
# 
# 
# 
# 
# ## Try a few other known cases
# ##
# 
# 
# 
# 








# c <- 137
# Etrans <- 3000 #eV
# Etrans_J <- Etrans_SI_keV *e
# Etrans_AU <- Etrans_SI_J/Eh
# ##other way
# E_AU <- 110
# E_J <- E_AU * Eh
# E_eV <- E_J / e
# ## ... works
# p <- Etrans_AU / c
# Spdb_AU <- compute_S(5,0.5,6,3400,Ar)





# The p bit
# m <- 1
# E <- p/c in atomic units 
# # stay in atomic units then one big conversion at the end


    
### ... ideally
#A0N0_S_pdb <- compute_S(A0N0_S_pdb)

##compare these to palffys theoretical estimates in her thesis using a well reasoned fraction

  
  #filter these observations where Ed - Ei = Egam => keep IC    .... A0N0_S_GAM$Etrans <- round(A0N0_S_GAM$Ed - A0N0_S_GAM$Ei, 2)
  
  # ### Current semi shit way of doing it
  # for(i in 1:length(A0N0_S_GAM$Etrans[i])){
  #   #if the transition matches assign A0N0_S$IC[i] <- CORRECT ND_Gamma$CC[i]
  #   
  #   Etrans <- round((A0N0_S$Ed[i] - A0N0_S$Ei[i]), 2)
#   if( Etrans == ND_Gamma$E[i] ) {
#     A0N0_S$IC[i] <- A0N0_S$CC[i]
#   }
#   else {
#     A0N0_S$IC[i] <- NA
#   }
#   
# }
# 
# A0N0_S <- filter(!is.na(A0N0_S$IC))

  


# #### Returns a value TRUE of FALSE if the IC coeff is known for the NEEC transition
# gammaTrans <- function(Etrans, Egam) {
#   if( Etrans == Egam ){
#     return(TRUE)  
#   }
#   else {
#     return(FALSE)
#   }
# }
# 
# A0N0_S_GAM$IC_NEEC <- map2(A0N0_S_GAM$Etrans, A0N0_S_GAM$Egam, gammaTrans)
# 






#### Taken from AnNi, finds and produces a vector of indexes
# for(i in 1:length(ND_AnNi$E)){
#   #store current isomer to be looked at
#   #progress(i)
#   Ax <- ND_AnNi$AX[i]
#   
#   #vector of isomer indexes
#   index <- which(isomer$AX==Ax)
# 
# }


###Assign IC column 
##Need to do another join here
#A0N0_S$IC <- filter(ND_Gamma, AX = ax)
# # try for loop
# for( i in 1:length(A0N0_S$IC) ) {
#   #if the transition matches assign A0N0_S$IC[i] <- CORRECT ND_Gamma$CC[i]
#   
#   #could use arr <-  which(ND_Gamma$AX == "12B")
#   #outputs : arr <- c(2 3 4)
#   
#   
#   #...shit way of doing it, nested for loops
#   #round(x, digits = 0)
#   for( j in 1:length(ND_Gamma)){
#     Etrans <- round((A0N0_S$Ed[i] - A0N0_S$Ei[i]), 2)
#     if( Etrans == ND_Gamma$E[j] ) {
#       A0N0_S$IC[i] <- ND_Gamma$CC[j]  
#     }
#   }
# }

  #ax <- levels(as.factor(ND_Gamma$AX))


 
