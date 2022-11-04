library(tidyverse)
library(plotly)
#for the dfactorial function
library(phangorn)
library(stringr)
library(plotly)
library(XML)
library(processx)
library(pbapply)

#### conv for addition of Z ####
Conv <- read.csv("Dependencies/Z_Conversion.csv")
conv <- data.frame(Z=Conv[1], Symb=Conv[2], Name=Conv[3])
conv$Symb <- toupper(conv$Symb)

source("Wins.R")

# ###Calculated Data
wap_N0 <- read_csv("App_Data/All_RatesCalcd_Viva.csv", col_types = cols( Ji = col_character(), ICC = col_double(), ICC_tot = col_double(), S = col_double(), S_tot = col_double(), N_neec_pw_Stot = col_double(), N_neec_pw = col_double()   ))
# #wap_Ni <- read_csv("Outputs/Wins_all_plasma_Ni.csv", col_types = cols( Ji = col_character(), ICC = col_double(), ICC_tot = col_double(), S = col_double(), S_tot = col_double(), N_neec_pw_Stot = col_double(), N_neec_pw = col_double()   ))
# wap <- wap_N0

### START HERE PRESENTLY
case_file <-"Microscopic/Theory_PDB_GS.csv"

#source("PDB_AnN0.R")
#function declaration, no manipulation 
source("Microscopic/PDB_ResonanceStrength_Calculator.R")
source("Microscopic/WeisskopfCalculator.R")
source("Microscopic/ICC_Eval.R")
#source("PDB_Rates_Calculator.R")
T_PDB <- read_csv(case_file)

#isomer <- read_csv("Dependencies/isomer.csv")
Ebinds <- read.csv("Dependencies/Ebinds/NIST_Ebinds.csv", header=TRUE, fill=TRUE)

# ## Filtering for comparison
# T_PDB_57FE <- filter(T_PDB, AX == "57FE")
# wap_57FE <- filter(wap_N0, AX == "57FE", Occ == 1, Ef < 14.5)
# 
# wap_57FE <- mutate(wap_57FE, Sp_Theoretical = double(1), Residual_Sth = double(1))
# test <- wap_57FE

## Organise the column structure
## currently tau is the total mean lifetime
T_PDB <- select(T_PDB, AX, Z, El, A, Q, Vi, Vi_q0, Ee, Ji, Jd, Jat, NEECX, Occ, Div, Type, tau_SI, Gamma_eV, oneoveroneplusalpha, AroverGamma, B_wu, B_wu_2, Yneec_th_SI, Yneec_calc_SI, Ar_SI, Ar_AU, IC_P, IC_tot, pi2_p2, pi2_p2_th, Sth, Ref, Spdb_P, Spdb_tot, Ratio_P, Rneec )

### Handy Functions
#"a" is which observation in conv$Symb
whi <- function(a) {return(which(conv$Symb==a))}

# which_S <- function(T_PDB, inp_AX) {
#  return(filter( T_PDB , AX == inp_AX)$S )
#   
# }
# test$Sp_Theoretical <- which_S(T_PDB)
n_index <- read_csv("Dependencies/n_to_ICshell_Conversion.csv")
## New columns for increasing accuracy of Spdb_p
T_PDB <- mutate(T_PDB, Vratio = double(1), nmax = double(1), nh = double(1), S_p_search = double(1) )
#T_PDB <- filter(T_PDB, NEECX != TRUE)


## Compute scaled ICC value!
for(i in 1:length(T_PDB$AX)){
                                              #function(Type, Zin, Q, n, OUTPUT, MIX, MR, MIX_type){
    T_PDB$IC_P[i] <- Compute_ICC_singles(T_PDB$Type[i], T_PDB$Z[i], T_PDB$Q[i], (T_PDB$Occ[i]+ 1), "p", FALSE, 1, "M1+E2")
    #T_PDB$Vi_q0[i] <- T_PDB$Q[i] - (Compute_ICC_singles(T_PDB$Type[i], T_PDB$Z[i], T_PDB$Q[i], (T_PDB$Occ[i]+ 1), "Eic"))
    #try(T_PDB$Vi_q0[i] <- T_PDB$Q[i] - (Compute_ICC_singles(T_PDB$Type[i], T_PDB$Z[i], T_PDB$Q[i], (T_PDB$Occ[i]+ 1), "Eic", FALSE, 1, "M1+E2")))
}
  
  
#pi2_p2 column is in AU
for(i in 1:length(T_PDB$AX)){
  T_PDB$Ar_SI[i] <- compute_Ar(T_PDB$Type[i], T_PDB$A[i], T_PDB$Q[i], T_PDB$B_wu[i])
  T_PDB$Yneec_calc_SI[i] <- T_PDB$Ar_SI[i] * T_PDB$IC_P[i] *  ( (2*T_PDB$Jd[i] + 1)*(2* T_PDB$Jat[i] + 1) / (2*(2*T_PDB$Ji[i] + 1)) )
 
  Yn_AU <- T_PDB$Yneec_th_SI[i] * time_conv
  T_PDB$pi2_p2_th[i] <- (( T_PDB$Sth[i] / (((length_conv^2)*Eh_eV)*(1/barn) )  ) / Yn_AU) * (1 + T_PDB$IC_P[i])
  
  ## The continuum electron energy into AU
  Ee_J <- 1E3*T_PDB$Ee[i] * e
  Ee_AU <- Ee_J/Eh
  c <- 137.036
  p <- sqrt( ((Ee_AU + c^2)^2) - (c^4) ) / c
  T_PDB$pi2_p2[i] <- ( (2*(pi^2) ) / (p^2) )
  
  T_PDB$Vratio <- T_PDB$Vi / T_PDB$Vi_q0
  
  # T_PDB$Spdb_P[i] <- compute_S(T_PDB$Jd[i],T_PDB$Jat[i],T_PDB$Ji[i],T_PDB$IC_P[i], T_PDB$Ar_SI[i], T_PDB$Q[i], T_PDB$Vi[i])
  # T_PDB$Ratio_P[i] <- T_PDB$Spdb_P[i] / T_PDB$Sth[i]
  T_PDB$Branching_p[i] <- 1/( 1 + (T_PDB$IC_P[i]*T_PDB$Div[i]) )
  
  # 
  T_PDB$Spdb_P[i] <- compute_S(T_PDB$Jd[i],T_PDB$Jat[i],T_PDB$Ji[i],T_PDB$IC_P[i], T_PDB$Ar_SI[i], T_PDB$Q[i], T_PDB$Vi[i])
  if(T_PDB$Occ[i] == 0){
   T_PDB$Spdb_P[i] <- T_PDB$Spdb_P[i] * T_PDB$Vratio[i]
  }
  if(T_PDB$Occ[i] != 0){
     T_PDB$Spdb_P[i] <- T_PDB$Spdb_P[i] * T_PDB$Branching_p[i] * T_PDB$Vratio[i]
  }
  T_PDB$Ratio_P[i] <- T_PDB$Spdb_P[i] / T_PDB$Sth[i]
  

  

 }


#write.table(T_PDB, file=case_file, append=T, row.names=F, col.names=T,  sep=",")

# for(i in 1:length(T_PDB$Eratio)){
#   i <- 1
#   element <- filter(Ebinds, Z==whi(T_PDB$Z[i])) 
#   T_PDB$Eratio[i] <-  element$IE[length(element[,1])] / element$IE[1]
# }



# 
# ### Handy Functions
# #Load Z column into NNDC
# 
# for(i in 1:length(T_PDB$Z)){
#   T_PDB$El[i] <- paste(unlist(str_remove_all(T_PDB$AX[i], "\\d" )), collapse = "")
#   T_PDB$Z[i] <- whi(T_PDB$El[i])
# }


# ## Produce the E and J values from the raw database
# T_PDB$Ei <- double(1)
# T_PDB$Ji <- character(1)
# T_PDB$Ed <- double(1)
# T_PDB$Jd <- character(1)
# 
# # ## Assign Ji from GS database
# for(i in 1:length(T_PDB$Ji)){
# #Assign J0 by infalliby running through GD
#   print(i)
#   if(!is.na(GS$Ji[i])){
#     T_PDB$Ji[i] <- GS$Ji[ which(GS$AX == T_PDB$AX[i] ) ]
#   }
# }
# 
# ## OR Assign Ji from isomer database **** WILL NEED TO sort out multiple isomers in each nuclide
# for(i in 36:42){
#   # if(!is.na(isomer$J[i])){
#     #Assign Ji by infalliby running through isomer
#     T_PDB$Ji[i] <- isomer$J[ which(isomer$AX == T_PDB$AX[i] ) ]
#   #}
# }
 
 
 
 

## Assign Jd from ND databse
# d <-1
# T_PDB$E <- round(T_PDB$Ed, d)
# ND$E <- round(ND$E, d)
# for(i in 35:length(T_PDB$Jd)){
# # i <- 3
#    #if(!is.na(ND$J[i])){
#      T_PDB$Jd[i] <- ND$J[ which( (ND$AX == T_PDB$AX[i]) & (ND$E == T_PDB$Etrans[i])  ) ]
# 
#    #   #}
#  print(i)
# }
# #write.table(T_PDB, file="Theory_PDB.csv", append=T, row.names=F, col.names=T,  sep=",")


## Assigne Etrans for the isomeric cases Ed - Ei
#T_PDB$Etrans <- T_PDB$Ed - T_PDB$Ei





### Spdb_P ###
#### Raw data ####
## Make the Jpi strings into just J doubles
# for(i in 1:length(T_PDB$Jd)) {
#   T_PDB$Jd[i] <- paste(unlist(str_extract_all(T_PDB$Jd[i], "[0-9/]")), collapse = "")
#   T_PDB$Jd[i] <- eval(parse(text=T_PDB$Jd[i]))  
#   
#   T_PDB$Ji[i] <- paste(unlist(str_extract_all(T_PDB$Ji[i], "[0-9/]")), collapse = "")
#   T_PDB$Ji[i] <- eval(parse(text=T_PDB$Ji[i])) 
# }
# T_PDB$Jd <- as.numeric(as.character(T_PDB$Jd))
# T_PDB$Ji <- as.numeric(as.character(T_PDB$Ji))
#   


## Make IC_P readable and numeric
  # for(i in 1:length(T_PDB$IC_P)){
  #   T_PDB$IC_P[i] <- eval(parse(text=T_PDB$IC_P[i])) 
  # }
  # T_PDB$IC_P <- as.numeric(as.character(T_PDB$IC_P))


#T_PDB <- T_PDB[,c(1:7,11,9,10,8,12:26)]
  
# ## IC_tot
# for(i in 1:length(T_PDB$Spdb_tot)) {
#   T_PDB$Spdb_tot[i] <- compute_S(T_PDB$Jd[i], T_PDB$Jat[i], T_PDB$Ji[i], T_PDB$IC_tot[i], T_PDB$Ar_SI[i], T_PDB$Etrans[i], T_PDB$Ebind[i])
# }
# 
# T_PDB$Ratio_tot <- T_PDB$Spdb_tot / T_PDB$Sth
#rm(T_PDB)
#T_PDB <- read_csv("Theory_PDB.csv")




### Check theoretical Yneecs first
# T_PDB$Sth <- double(1)
# T_PDB <- T_PDB[,c(1:20,30,21:29)]

# T_PDB$Gamma_AU <- T_PDB$Gamma / Eh_eV
# T_PDB <- T_PDB[,c(1:15,31,16:30)]
# 
# for(i in 1:length(T_PDB$Sth_Check)) {
#   print(i)
#   T_PDB$Sth[i] <- compute_S(T_PDB$Jd[i], T_PDB$Jat[i], T_PDB$Ji[i], T_PDB$IC_P[i], T_PDB$Yn[i], T_PDB$Etrans[i], T_PDB$Ebind[i])
# }
#T_PDB$Ratio_P <- T_PDB$Spdb_P / T_PDB$Sth  

#### Ar_AU/Gamma_Au ####
#T_PDB$oneoveroneplusalpha <- double(1)
T_PDB$oneoveroneplusalpha <- 1/(1+T_PDB$IC_P)
# T_PDB <- T_PDB[,c(1:16,32,17:31)]
# 
#T_PDB$ArAUoverGammaAU <- double(1)
T_PDB$ArAUoverGammaAU <- T_PDB$Ar_AU / T_PDB$Gamma_AU
# T_PDB <- T_PDB[,c(1:17,33,18:32)]
# 
# #### Must convert from Bdowns to Bups!! ###
# T_PDB <- T_PDB %>% mutate(B_wu_up = (((2*Jd) + 1) / ((2*Ji) + 1)) * B_wu)
# T_PDB <- T_PDB[,c(1:19,34,20:33)]

# T_PDB$Ratio_SanityCheck <- T_PDB$ArAUoverGammaAU / T_PDB$oneoveroneplusalpha
# T_PDB <- T_PDB[,c(1:18,36,19:35)]


###WORKING UP TO HERE
#write.table(T_PDB, file=case_file, append=T, row.names=F, col.names=T,  sep=",")
# unlink(case_file)
# write.table(T_PDB, file=case_file, append=T, row.names=F, col.names=T,  sep=",")
# T_PDB <- read_csv(case_file)


# ## Compute Ar with B_wu_up
# T_PDB$Ar_SI <- 1
# for(i in 1:length(T_PDB$Ar_SI)){
#   T_PDB$Ar_SI[i] <- compute_Ar(T_PDB$Type[i], T_PDB$A[i], T_PDB$Etrans[i], T_PDB$B_wu[i]) 
# }
# 
# T_PDB$Ar_AU <- T_PDB$Ar_SI * time_conv
# 
# 
# ##check the theoretical branching ratios against my calcultion of (1/1+alpha)
# 
# 
# ## IC_P
# for(i in 1:length(T_PDB$Spdb_P)) {
#   print(i)
#   T_PDB$Spdb_P[i] <- compute_S(T_PDB$Jd[i], T_PDB$Jat[i], T_PDB$Ji[i], T_PDB$IC_P[i], T_PDB$Ar_SI[i], T_PDB$Etrans[i], T_PDB$Ebind[i])
# }
# 
# T_PDB$Ratio_P <- T_PDB$Spdb_P / T_PDB$Sth
# 
# 
# 
# ## COmpute EBIT rate
# 
# ## compute_S(T_PDB$Jd[i], T_PDB$Jat[i], T_PDB$Ji[i], T_PDB$IC_K[i], T_PDB$Ar_SI[i], (T_PDB$Etrans[i] * 1E3)) 
# 
# #find the right B to use in the Ar calc
# 
# ##Finish proper Ar rate calculator
# 
# ##Fill in missing IC values with BRICC
# 
# ##Fill in the S PDB tot and K columns
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
# 
# 
# # ## testing Pallfy Ar computation 84Rb
# # En_MeV <- 3.05E-3
# # BW <- 0.08
# # B_SI <- 1.790 * BW
# # Ar_SI <- 1.779E13*(En_MeV^3)*B_SI
# # Ar_SI #[s^-1]
# #   
# # 
# # #Ar_AU_conv == Ar_AU_calc
# # # time_AU <- seconds * time_conv
# # # rate_AU <- rate_SI / time_conv
# #  #Ar_SI <- Ar_AU * time_conv
# # B_AU <- B_SI * 0.25 * (me.mp)^2 #* 1.85480201566E23
# # En_AU <- (En_MeV * 1E6)*(e/Eh)
# # Ar_AU_calc <- compute_Ar(1, B_AU, En_AU  )
# # Ar_SI_calc <- Ar_AU_calc * time_conv ## WRONG UPON EXAMINATION
# # Ar_SI_calc
# # 
# # ## OK go for the table read for now
# 
# 
# 






