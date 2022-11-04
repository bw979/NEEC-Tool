source("Wins.R")

##########################
######## PLASMA ##########
##########################
FLYCHK_ne <- c(1E12, 1E13, 1E14, 1E15, 1E16, 1E17, 1E18, 1E19, 1E20, 1E21, 1E22, 1E23, 1E24)
FLYCHK_Te <- c(0.5, 1, 1.5, 2, 5, 7, 10, 15, 23, 32, 52, 74, 100, 165, 235, 310, 390, 475, 655, 845, 1000, 1441, 1925, 2454, 3030, 3655, 4331, 5060, 5844, 6685, 7585, 8546, 10000, 20000, 50000, 100000)
##For the subshell column
ICC_Multipliers <- read_csv("Dependencies/ICC_Multipliers.csv") 
#ICC_Multipliers <-  arrange(ICC_Multipliers, Z, Occ)
atomic_colours <- read_csv("Dependencies/Atomic_Colours.csv")

# source("Macroscopic/PDB_Rates_Plasma.R")
# source("Microscopic/ICC_Eval.R")

# for(n in 1:104){
# 
#   Name_AnN0_Rates <- sprintf("Outputs/MASTER/A%dN0_Plasma_Rates.csv", n)
#   Name_AnNi_Rates <- sprintf("Outputs/MASTER/A%dNi_Plasma_Rates.csv", n)
#     
# unlink(Name_AnN0_Rates)
# unlink(Name_AnNi_Rates)
# 
# }



##### THIS FILE FOR FINAL RUN ... COLLECT AND RENAME THE WRANGLED n files in MASTER/All/ ####

##### GROUND STATE #####
# Read and merge all An data sets and all Ni datasets
# Merge all data AnN0 together
# AnN0_all <- list.files(path = "Outputs/MASTER/",     # Identify all csv files in folder
#                        pattern = "(A\\d+N0_Plasma_Rates.csv)", full.names = TRUE) %>%
#   lapply(read_csv) %>%                                            # Store all files in list
#   bind_rows             # Combine data sets into one data set
# 
# if(!file.exists("Outputs/MASTER/All/AnN0_all.csv")){
# write.table(AnN0_all, file="Outputs/MASTER/All/AnN0_all.csv", append=T, row.names=F, col.names=T,  sep=",")
# }


# AnN0_all <- read_csv("Outputs/MASTER/All/AnN0_all.csv")
# AnN0_all <- mutate(AnN0_all, Isomer_Egam = NA, Isomer_Decay_Type = NA)
# Wins_all_location_save <- "Outputs/MASTER/AnN0_all_RatesCalcd.csv"
# #AnNi_all <- read_csv("Outputs/MASTER/All/AnNi_all.csv")
# 
#AnNi_all<-read_csv("Outputs/MASTER/All/Wins_all_plasma_Ni.csv")
# 
# # ###### ISOMERIC ####
# # AnNi_all <- list.files(path = "Outputs/MASTER/",     # Identify all csv files in folder
# #                        pattern = "(A\\d+Ni_Plasma_Rates.csv)", full.names = TRUE) %>%
# #   lapply(read_csv) %>%                                            # Store all files in list
# #   bind_rows                                                    # Combine data sets into one data set
# # if(!file.exists("Outputs/MASTER/All/AnNi_all.csv")){
# #   write.table(AnNi_all, file="Outputs/MASTER/All/AnNi_all.csv", append=T, row.names=F, col.names=T,  sep=",")
# # }
# 
# 
# ########Filter for distinct: priority order AX, Ef, Thf %>% !is.na(Thf) ... 1E3 HORRAY
# AnN0_all_dist <- arrange(AnN0_all, M, AX, desc(Occ)) %>%
#   distinct(AX,Ef, Occ, .keep_all = T) %>%
#   filter(!is.na(B))


# AnNi_all <- read_csv("Outputs/MASTER/All/AnNi_all.csv")
# Wins_all_location_save <- "Outputs/MASTER/AnNi_all_RatesCalcd.csv"
# 
# ########Filter for distinct: priority order AX, Ef, Thf %>% !is.na(Thf) ... 1E3 HORRAY
# AnNi_all_dist <- arrange(AnNi_all, M, AX, desc(Occ)) %>%
#   distinct(AX, Ei, Ef, Occ, .keep_all = T) %>%
#   filter(!is.na(B), !is.na(Thf_numeric))
# 
# #Wap <- bind_rows(AnN0_all_dist, AnNi_all_dist)
# #Wap <- AnN0_all_dist
# Wap <- AnNi_all_dist
# 
# ## Add mixing columns and ICC-tot
# Wap <- Wap %>%
#   mutate(Jf_double=double(1), Ji_double=double(1), MIX_type=character(1), ICC_tot=double(1), Stot=double(1),
#          Gamma_Thalf_eV=double(1), Gamma_eV=double(1), ICC_sum=double(1), Npw_Pl_S = double(1), S_sum=double(1), Npw_SEBIT=double(1) )
# 
# #### FINAL DATABASE STRUCTURE
# Wap <- select(Wap, AX, M, El, Z, N, Ei, Ef, Ji, Ji_double, Jf, Jf_double, Type, Thi, Thi_numeric, Thf, Thf_numeric,
#                Occ, CS, At_GS_Config, subshell, Jat, Ebind, Q, Egam, Efeed, FL, Btype, B, MR, MIX_type, B2type, B2, RI, Ar, Gamma_Thalf_eV, Gamma_eV,
#               ICC, ICC_tot, ICC_sum, S, Stot, S_sum, Ee, Ebeam, -Rate, Npw_Pl_S, Npw_SEBIT)

# Wap <- Wap %>% filter(M >= m1 & M <= m2)

# ## For J columns
# # Wap$Jf_double <- double(1)
# # Wap$Ji_double <- double(1)


#### 178HF example #########
#Wap <- read_csv("Outputs/MASTER/All/Wap_All.csv")
########EBIT FILTER
#Wap_dist <- arrange(Wap, M, AX, Occ) %>%
#  filter(AX=="178HF", !is.na(B)) 
#Wap <- Wap_dist


#### USE THIS BIT IF WE HAVE ALREADY HAD A RUN OF RATES
#Wap <- read_csv("Outputs/MASTER/All/All_RatesCalcd.csv")
#Wap <- read_csv("Outputs/MASTER/All/All_RatesCalcd_Viva.csv")
Wap <- read_csv("Mo_Additions_Rates.csv")
#Wap2 <- Wap[,c(1:23,43, 44, 24:42,45:46)] 
#write.table(Wap2, file="Outputs/MASTER/All/All_RatesCalcd.csv", append=T, row.names=F, col.names=T,  sep=",")

#Wap_tot <- distinct(Wap, Ef, .keep_all = T)

#Wap <- select(Wap, -Rneec_Pl, Rneec_SEBIT)
Wap$Q <- Wap$Ef - Wap$Ei
Wap$Ee <- Wap$Q - Wap$Ebind
Wap <- filter(Wap, Ee>0)
#Wap$Ebeam <- mp.me *


##### PRESENT... WANT TO calculate ICC's and do some initial rate calcs
for(i in 1:nrow(Wap) ){
#for(i in 1:100 ){
  print(i)
  
   
  # # ### some extra wrangling to add subshell column (xray notation)
  # # #Wap$subshell[i] <- filter(ICC_Multipliers, Z==Wap$Z[i], Occ==Wap$Occ[i])$subshell
  # # Wap$shell[i] <- str_remove_all(Wap$subshell[i], "\\d")
  # # Wap$shell_colour[i] <- atomic_colours$colour2[which(atomic_colours$subshell==Wap$subshell[i])[1]] 
  # # 
  # 
  ### Create the mixing type column eg. "M1+E2"
  if( !is.na(Wap$Btype[i]) && !is.na(Wap$B2type[i]) && Wap$Btype[i] =="BM1W" && Wap$B2type[i]=="BE2W=" && Wap$MR[i]!=0 && !is.na(Wap$MR[i])  ){
  #if( !is.na(Wap$Btype[i])  && Wap$Btype[i] =="BM1W" && Wap$B2type[i]=="BE2W=" && Wap$MR[i]!=0 && !is.na(Wap$MR[i])  ){
    Wap$MIX_type[i] <- "M1+E2"
  } else {
    Wap$MIX_type[i] <- ""
  }

  ###### CONVERT Jpi to just numeric J values
  if(Wap$Ji[i]=="+") Wap$Ji_double[i] = NA
  if(Wap$Jf[i]=="+") Wap$Jf_double[i] = NA
  if(Wap$Ji[i]=="-") Wap$Ji_double[i] = NA
  if(Wap$Jf[i]=="-") Wap$Jf_double[i] = NA

  if(Wap$Ji[i]!="+" && Wap$Jf[i]!="+" && Wap$Ji[i]!="-" && Wap$Jf[i]!="-"){

    # Make the Jpi strings into just J doubles
    Wap$Jf_double[i] <- eval(parse(text=paste(unlist(str_extract_all(Wap$Jf[i], "[0-9/]")), collapse = "")))
    Wap$Ji_double[i] <- eval(parse(text= paste(unlist(str_extract_all(Wap$Ji[i], "[0-9/]")), collapse = "")  ))
    # Wap$Jf_double[i] <- as.double(as.character(Wap$Jf_double[i]))
    # Wap$Ji_double[i] <- as.double(as.character(Wap$Ji_double[i]))
  }

  ##### CALCULATE Ar, ICC and Stot
  if(  !is.na(Wap$B[i]) && !is.na(Wap$Ji_double[i]) && !is.na(Wap$Jf_double[i]) && !is.na(Wap$Type[i]) && !is.na(Wap$ICC[i]) && Wap$Type[i] != "E0" && Wap$Type[i] != "Large_L") {
    ######Compute_Ar... reuturns in SI [s^-1] ... also trying ICC_tot
    #if its mixed
    if(  (Wap$MR[i]!=0) && (!is.na(Wap$MR[i])) ){
      Wap$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i], TRUE, Wap$MR[i] )
      ## If its "M1 + E2" mixed
      if( !is.na(Wap$MIX_type[i]) && Wap$MIX_type[i] == "M1+E2"  ){
        try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot", TRUE,  Wap$MR[i], Wap$MIX_type[i])  )
        try(Wap$ICC[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", TRUE, Wap$MR[i], Wap$MIX_type[i])  )
        try(Wap$ICC_sum[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "SUMp", TRUE, Wap$MR[i], Wap$MIX_type[i])  )
      }
    } else {
      Wap$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i], FALSE, Wap$MR[i] )
      try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot", FALSE,  0, ""))
      try(Wap$ICC[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", FALSE,  0, ""))
      # try(Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot"))
      # try(Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot"))
      try(Wap$ICC_sum[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "SUMp", FALSE, 0, "")  )
    }
  }

  
  #   # #Compute_S Jd, Jat,Ji,ICC, Ar, Etrans_keV, Ebind_keV
  #   # Wap$S[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
  #   # Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
  # 
  # 
  # ### Total depletion level width
  Wap$Gamma_Thalf_eV[i] <- (hbar/e) * (log(2) / Wap$Thf_numeric[i])
  Wap$Gamma_eV[i] <-  (hbar/e) * (Wap$Ar[i]*(1+Wap$ICC_sum[i]))
  Wap$Gamma_Scale_Factor[i] <- 1/(1+Wap$ICC_sum[i])

  #### RATES ####

  j <- which(facility$Facility == "ELI-NP")   #...ELI NP
  
    ### COMPUTE THE RESONANCE STRENGTH

    Wap$Stot[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], mean(Wap$Ebind))*Wap$Z[i]
    Wap$S[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i]) *Wap$Gamma_Scale_Factor[i]
    Wap$S_sum[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_sum[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])


  
  #if(Z<4 ||)Z>79) {Wap$Npw_Pl_Stot_Eli[i]}
 
  try(Wap$Rate_Pl[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], 100E-12, "Te", (1000*facility$Te_Beg_keV[j]), "Rate", "p"))
  try(Wap$Rate_SEBIT[i] <- compute_EBIT_RATE(Wap$S[i], Wap$Ee[i], Wap$Gamma_eV[i], "Rate") )
  try(Wap$Npw_Pl_S[i] <- compute_Plasma_Rate(Wap[i,], facility[j,] , 100E-12, "Te", (1000*facility$Te_Beg_keV[j]), "N_neec_pw", "p"))
  try(Wap$Npw_SEBIT[i] <- compute_EBIT_RATE(Wap$S[i], Wap$Ee[i], Wap$Gamma_eV[i], "N_neec_pw") )
  
}

#Wins_all_location_save <- sprintf("Outputs/MASTER/Chunk/Wins_all_plasma_%003.f_%003.f.csv", m1, m2)
Wins_all_location_save <- "Mo_Additions_Rates2.csv"

if(!file.exists(Wins_all_location_save)){
  write.table(Wap, file=Wins_all_location_save, append=T, row.names=F, col.names=T,  sep=",")
}





# ######### NOW RUN IC, Ar, and S calculator..... including mixing if available
# #for(i in 1:nrow(Wap) ){
# for(i in 1:500 ){
#  print(i)
#   
#   # #### Create the mixing type column eg. "M1+E2"
#   # if( !is.na(Wap$Btype[i]) && !is.na(Wap$B2type[i]) && Wap$Btype[i] =="BM1W" && Wap$B2type[i]=="BE2W=" && Wap$MR[i]!=0 && !is.na(Wap$MR[i])  ){
#   #   Wap$MIX_type[i] <- "M1+E2"
#   # } else {
#   #   Wap$MIX_type[i] <- ""
#   # }
#   
#   ###### CONVERT Jpi to just numeric J values
#   if(Wap$Ji[i]=="+") Wap$Ji_double[i] = NA
#   if(Wap$Jf[i]=="+") Wap$Jf_double[i] = NA
#   if(Wap$Ji[i]=="-") Wap$Ji_double[i] = NA
#   if(Wap$Jf[i]=="-") Wap$Jf_double[i] = NA
#   
#   if(Wap$Ji[i]!="+" && Wap$Jf[i]!="+" && Wap$Ji[i]!="-" && Wap$Jf[i]!="-"){
#     
#     # Make the Jpi strings into just J doubles
#     Wap$Jf_double[i] <- eval(parse(text=paste(unlist(str_extract_all(Wap$Jf[i], "[0-9/]")), collapse = "")))
#     Wap$Ji_double[i] <- eval(parse(text= paste(unlist(str_extract_all(Wap$Ji[i], "[0-9/]")), collapse = "")  ))
#     # Wap$Jf_double[i] <- as.double(as.character(Wap$Jf_double[i]))
#     # Wap$Ji_double[i] <- as.double(as.character(Wap$Ji_double[i]))
#   }
#   
#   ##### CALCULATE Ar, ICC and Stot
#   if(  !is.na(Wap$B[i]) && !is.na(Wap$Ji_double[i]) && !is.na(Wap$Jf_double[i]) && !is.na(Wap$Type[i]) && Wap$Type[i] != "E0" && Wap$Type[i] != "Large_L") {
#     ######Compute_Ar... reutrns in SI [s^-1] ... also trying ICC_tot
#     #if its mixed 
#     if(  (Wap$MR[i]!=0) ){
#       Wap$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i], TRUE, Wap$MR[i] )
#       ## If its "M1 + E2" mixed
#       try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot", FALSE,  0, ""))
#       try(Wap$ICC[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", FALSE,  0, ""))
#     } else {   
#     #if((!is.na(Wap$MR[i])) &&) {
#       Wap$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i], FALSE, Wap$MR[i] )
#       try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot", FALSE,  0, ""))
#       try(Wap$ICC[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", FALSE,  0, ""))
#       #try(Wap$ICC_sum[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "SUMp", FALSE, 0, "")  )
#     }
#     
#       if(  Wap$MIX_type[i] == "M1+E2"  ){
#         try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot", TRUE,  Wap$MR[i], "M1+E2") )
#         try(Wap$ICC[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", TRUE, Wap$MR[i], "M1+E2")  )
#         #try(Wap$ICC_sum[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "SUMp", TRUE, Wap$MR[i], "M1+E2")  )
#       }
#     
#     
#     ### COMPUTE THE RESONANCE STRENGTH
#     Wap$Stot[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#     Wap$S[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#     #Wap$S_sum[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_SUM[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#     
#     # #Compute_S Jd, Jat,Ji,ICC, Ar, Etrans_keV, Ebind_keV
#     # Wap$S[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#     # Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#   }
#   
#   ### Total depletion level width
#   Wap$Gamma_Thalf_eV[i] <- (hbar/e) * (log(2) / Wap$Thf_numeric[i])
#   #Wap$Gamma_eV[i] <-  (hbar/e) * (Wap$Ar[i]*(1+Wap$ICC_SUM[i]))
#   Wap$Gamma_eV[i] <-  (hbar/e) * (Wap$Ar[i]*(1+Wap$ICC_tot[i]))
#   
#   #### RATES ####
#   j<-21
#   #if(Z<4 ||)Z>79) {Wap$Npw_Pl_Stot_Eli[i]}
#   try(Wap$Npw_Pl_S[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", (1000*facility$Te_Beg_keV[j]), "N_neec_pw", "p"))
# }
# 
# # if(!file.exists("Outputs/MASTER/All/WAP.csv")){
# #   write.table(Wap, file="Outputs/MASTER/All/WAP.csv", append=T, row.names=F, col.names=T,  sep=",")
# # }

###For calculating the NEEC width
# Gamma_tot_eV <- log(2) * (hbar / Wap$Thf_numeric[i]) *(1/e)
# Gamma_d_eV <- hbar * (Wap$Ar[i]*(1+Wap$ICC_tot[i])) * (1/e)

######### Now Run NEEC rates for ALL scenarios
# #Rneec_PF, Rneec_Ch_Argonne, Rneec_Ch_Triumf, Rneec_Isolde, Rneec_S_EBIT, Rneec_GSI


# Wap <- mutate(Wap, Npw_PL_Stot=double(1)) %>%
#   filter(Z<80)
#Wins_all, facility, INPUT_Te, Te_input, OUTPUT, OUTPUT2
# ## Choose facility
# j<- 16
# for(i in 1:nrow(Wap)){
#   print(i)
#   #i<-1000
#   Wap$Npw_PL_Stot[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", (1000*facility$Te_Beg_keV[j]), "N_neec_pw", "tot")
# }







######################## ATOMIC SPECTRUM ANALYSIS  ##################
# Wap_n <- AnN0_all %>% distinct(AX, Ef, Occ, .keep_all = T) %>%
#   arrange(M) %>%
#   filter(!is.na(B))
  















# check <- read_csv("Outputs/MASTER/A3N0_Plasma_Rates.csv")
# 
# # nuclide <- "84RB"
# # 
# # 
# #### Calculated data filter nuclide of interest
# wap_N0 <- read_csv("Outputs/Wins_all_plasma_N0.csv", col_types = cols( Ji = col_character(), ICC = col_double(), ICC_tot = col_double(), S = col_double(), S_tot = col_double(), N_neec_pw_Stot = col_double(), N_neec_pw = col_double()   ))  
#  # filter(AX == nuclide)
# #wap_Ni <- read_csv("Outputs/Wins_all_plasma_Ni.csv", col_types = cols( Ji = col_character(), ICC = col_double(), ICC_tot = col_double(), S = col_double(), S_tot = col_double(), N_neec_pw_Stot = col_double(), N_neec_pw = col_double()   )) 
#  # filter(AX == nuclide)
# # 
# Wap <- wap_N0
# # Wap <- select(Wap, -N_neec_pw_Stot, -N_neec_pw, -Rate)
# # 
# # write.table(Wap, file="Outputs/Wins_all_plasma_N0.csv", append=T, row.names=F, col.names=T,  sep=",")
# # 
# #Wap <- rbind(wap_N0, wap_Ni)
# # #Wap <- mutate(Wap, S_tot = double(1))
# # # CHUNK ND
# # 
# # 
# args <- commandArgs(TRUE)
# m1 <- as.numeric(args[1])
# m2 <- as.numeric(args[2])
# Wap <- Wap %>% filter(M >= m1 & M <= m2)
# # # 
# # #### Distinct levels ####
# # Dlevels <- distinct(Wap, Ef, Type, .keep_all = TRUE) %>%
# #   filter(!is.na(Thf_numeric)) %>%
# #   distinct(Thf_numeric, .keep_all = TRUE)
# # 
# # 
# # 
# # #### ADDING IN THE MIXING RATIO  #####
# # #-Round data for comparison
# # 
# # #- Add B2 and delta column
# # Wap <- mutate(Wap, B2 = double(1), delta = double(1))
# # 
# # #-Assign B2 and MR
# # ## Rount to 3sf
# # signif <- 3
# # Dlevels$Q <- signif(Dlevels$Q, signif)
# # Dlevels <- unique(Dlevels) %>% filter(!is.na(Jf), !is.na(Ji))
# # 
# # Gamma_Bvals <- filter(Gamma_Bvals, AX == nuclide)
# # Gamma_Bvals$Egam <- signif(Gamma_Bvals$Egam, signif)
# # 
# # # do the grab and assign the new variables
# # Gamma_Bvals <- mutate(Gamma_Bvals, B2 )
# # 
# # ## Regain precision in Q
# # Wap$Q <- Wap$Ee + Wap$Ebind
# # ####   #####
# 
# 
# ############# CSD  #############
# #- Add a CSD fraction for each of the Occ, At all the densities and electron temperatures another 468 columns
# 
# 
# 
# # - Wrangle the data structure here at least by tonight
# 
# # m1<-56;m2<-58
# # Wap <- Wap %>% filter(M >= m1 & M <= m2)
# 
# ############       #############
# 
# 
# 
# # S column means S_p
# n_index <- read_csv("Dependencies/n_to_ICshell_Conversion.csv")
# 
# ## ELI-NP j=16
# j<-16
# #registerDoParallel(15) 
# #foreach (i=1:length(Wap$ICC_tot)) %dopar% {
# Wap <- mutate(Wap, WE = FALSE, Gamma_neec_eV = double(1))
# 
# 
# for(i in 1:length(Wap$Ar)) {
#   print(i)
#   ## Add B=1 if NA and label with WE =TRUE
#    if( is.na(Wap$B[i]) ){
#      Wap$B[i] <- 1
#      Wap$WE[i] <- TRUE
#    }
#   
#   ## Recalculate Ar's and S's
#   if(  !is.na(Wap$B[i]) && !is.na(Wap$Ji_double[i]) && !is.na(Wap$Jf_double[i]) && !is.na(Wap$ICC[i]) && Wap$Type[i] != "E0" && Wap$Type[i] != "Large_L") {
#     ######Compute_Ar
#     Wap$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i])
#     #Compute_S Jd, Jat,Ji,ICC, Ar, Etrans_keV, Ebind_keV
#     Wap$S[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#     Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#   }
#   
#   Gamma_neec_J <-  hbar * (Wap$Ar[i]*(1+Wap$ICC_tot[i]))
#   Wap$Gamma_neec_eV[i] <- Gamma_neec_J / e
#   
#   if(!is.na(Wap$S[i])){
#   ##PLASMA RATE AND N_neec_pw
#   j <- 16  
#   Wap$N_neec_pw_Plasma_Stot_Eli[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", 1000*facility$Te_Beg_keV[j], "N_neec_pw", "tot")
#   Wap$N_neec_pw_Plasma_Eli[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", 1000*facility$Te_Beg_keV[j], "N_neec_pw", "p")
#   Wap$R_neec_Plasma_Eli[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", 1000*facility$Te_Beg_keV[j], "Rate", "p")
#   j <- 17  
#   Wap$N_neec_pw_Plasma_Stot_NdG[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", 1000*facility$Te_Beg_keV[j], "N_neec_pw", "tot")
#   Wap$N_neec_pw_Plasma_NdG[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", 1000*facility$Te_Beg_keV[j], "N_neec_pw", "p")
#   Wap$R_neec_Plasma_NdG[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Te", 1000*facility$Te_Beg_keV[j], "Rate", "p")
#   
#   }
#   #Wins_all_plasma$shell[i] <- n_index$Shell[Wins_all_plasma$Occ[i]]
#   ## Making S_n ... more accurate
#   
#   ## Compute SEBIT Rate
#   #compute_EBIT_RATE <- function(S, Ee, V, nion, nelec, s, Gamma_neec_eV, OUTPUT) {
#   if(Wap$Ee[i] < 260){
#     Wap$Rneec_SEBIT[i] <- compute_EBIT_RATE(Wap$S[i], Wap$Ee[i], V, nion, nelec, s, Wap$Gamma_neec_eV[i], "rate" )
#     Wap$N_neec_pw_SEBIT[i] <- compute_EBIT_RATE(Wap$S[i], Wap$Ee[i], V, nion, nelec, s, Wap$Gamma_neec_eV[i], "N_neec_pw" )
#   } else {
#     Wap$Rneec_SEBIT[i] <- NA
#     Wap$N_neec_pw_SEBIT[i] <- NA
#   }
# }
# 
# #Wins_all_location_save <- sprintf("Outputs/Chunk_N0/Wins_all_plasma_%003.f_%003.f.csv", m1, m2)
# #Wins_all_location_save <- sprintf("Outputs/Chunk_Ni/Wins_all_plasma_%003.f_%003.f.csv", m1, m2)
# write.table(Wap, file=Wins_all_location_save, append=T, row.names=F, col.names=T,  sep=",")
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
# # ## Rount to 2sf
# # figs <- 2
# # Gamma_Bvals$Egam <- signif(Gamma_Bvals$Egam, figs)
# # #Gamma_Bvals$Egam <- ceiling(Gamma_Bvals$Egam)
# # 
# # ### Add the 93Mo Bvalue with the new Gamma additions
# # for(i in 1:length(Wins_all$AX)){
# #   print("Fill in Known B values")
# #   print(i)
# #   ###Sort out B values
# #   if(Wins_all$AX[i] == "93MO"){
# #     G <- filter(Gamma_Bvals, AX == Wins_all$AX[i], Egam == Wins_all$Q[i], !is.na(Btype))
# #     if(length(G$AX)==0) Wins_all$B[i] = NA
# #     else Wins_all$B[i] <- G$B
# #   }
# # }
# # 
# # Wins_all_Bknown <- filter(Wins_all, !is.na(B), B!=0)
# # length(Wins_all_Bknown$B)
# 
# # Wins_all_location <- "Outputs/Wins_all.csv"
# # unlink(Wins_all_location)
# # if(!file.exists(Wins_all_location)){
# #   write.table(Wins_all, file=Wins_all_location, append=T, row.names=F, col.names=T,  sep=",")
# # }
# 
# ## Create iso_num column (zeroth isomer is ground state)
# Wins_all <- mutate(Wins_all, Unity_Plasma_Rate = double(1), N_neec_pw = double(1))
# 
# 
# #Wap <- filter( Wins_all, !is.na(IC), IC != 0, !is.na(B), B != 0  ) 
# ## set j which is increment through facility database
# # #Generic long pulse ND:Glass
# # j <- 1
# # ### Rates Calculation for loop
# # for(i in 1:length(Wins_all$S)){ 
# #   #i<-1
# #   
# #   print(i)
# #   ## EBIT RATE
# #   #AnN0$EBIT_Rate[i] <- compute_EBIT_RATE(AnN0$S[i], AnN0$Ee[i], V, nion, nelec, spread)
# #   
# #   #j <- 1
# #   #################### Needs J indexing properly
# #   #for( j in 1:length(facility$Irradiance_Wcm2um2) ){
# #     ## PLASMA RATE
# #     if( is.na(Wins_all$S[i]) || is.na(Wins_all$Ee[i]) || is.na(Wins_all$Thf_numeric[i]) ) { 
# #       Wins_all$N_neec_pw[i] = NA
# #       Wins_all$Unity_Plasma_Rate[i] = NA
# #     }
# #     else {
# #       #Wins_all$N_neec_pw[i] <- compute_Plasma_Rate(Wins_all[i,], facility[j,])
# #       Wins_all$N_neec_pw[i] <- compute_Plasma_Rate(Wins_all[i,], facility[j,], "N_neec_pw")
# #       Wins_all$Unity_Plasma_Rate[i] <- compute_Plasma_Rate(Wins_all[i,], facility[j,], "Rate")
# #     }
# #   #}
# #     
# #   ##### compute rate for each facility
# #   #Wins_all
# #   
# #   ## Charge state distribution rate
# # }
# 
# ###### CSD #######
# ##First Make a discrete probability distribution array
# # CS <- 10
# # m0 <- dnorm(10, mean = CS, sd = 3)
# # mm1 <- dnorm(9, mean = CS, sd = 3)
# # mm2 <- dnorm(8, mean = CS, sd = 3)
# # mm3 <- dnorm(7, mean = CS, sd = 3)
# ####### Just use CS = Z+ numbers for now
# 
# 
# 
# 
# # Wins_all_location <- "Outputs/Wins_all.csv"
# # if(!file.exists(Wins_all_location)){
# #   write.table(Wins_all, file=Wins_all_location, append=T, row.names=F, col.names=T,  sep=",")
# # }
# 
# ## Looking at the data trends
# #Wins_all2 <- filter( Wins_all, !is.na(Unity_Plasma_Rate), Unity_Plasma_Rate!=0 )
# #Wins_all2 <- filter( Wins_all, !is.na(Nexc_pw), Nexc_pw!=0 )
# #Wins_all_test <- filter( Wins_all, !is.na(Rate), Rate != 0) 
# #Wap <- filter( Wins_all, !is.na(N_neec_pw), N_neec_pw != 0) 
# #Wap <- filter( Wins_all, !is.na(Unity_Plasma_Rate), Unity_Plasma_Rate != 0) 
# Wap <- Wins_all
# #compute_Plasma_Rate(Wins_all_test[1,], facility[9,])
# 
# ##### Now add ICC_tot
# Wap <- mutate(Wap, ICC_tot = double(1)) %>%
#   filter(Z>11, !is.na(Type), Type != "E0")
# 
# # Wins_all_location <- "Outputs/Wap.csv"
# # unlink(Wins_all_location)
# # if(!file.exists(Wins_all_location)){
# #   write.table(Wap, file=Wins_all_location, append=T, row.names=F, col.names=T,  sep=",")
# # }
# Wins_all_location <- "Outputs/Wap.csv"
# Wap <- read_csv(Wins_all_location)
# Wap$Q <- Wap$Ee + Wap$Ebind
#Wap_test <- filter(Wap, AX == "57FE")
## Try the new ICC eval function

# ########## Subdivide the problem
# sublength <- length(Wap$ICC_tot) / 10 
# sub_round <- round(sublength)
#l_last <- length(Wap$ICC_tot)

# for(i in 378422:500000) {
#   Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot")
# }

# #Wap <- read_csv(Wins_all_location)
# for(i in 1:500000) {
#   Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot")
# }
# unlink(Wins_all_location)
# if(!file.exists(Wins_all_location)){
#   write.table(Wap, file=Wins_all_location, append=T, row.names=F, col.names=T,  sep=",")
# }
# #Wap <- read_csv(Wins_all_location)
# for(i in 500001:1000000) {
#   Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot")
# }
# unlink(Wins_all_location)
# if(!file.exists(Wins_all_location)){
#   write.table(Wap, file=Wins_all_location, append=T, row.names=F, col.names=T,  sep=",")
# }
#Wap <- read_csv(Wins_all_location)
##for(i in 1000001:1500000) {
#  v <- c(1000001:1500000)
#  #Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot")
#  Wap$ICC_tot[v] <- pmap(list(Wap$Type[v], Wap$Z[v], Wap$Q[v], "tot"), Compute_ICC_singles )
##}
#unlink(Wins_all_location)
#if(!file.exists(Wins_all_location)){
#  write.table(Wap, file=Wins_all_location, append=T, row.names=F, col.names=T,  sep=",")
##}

# ###########################################################################
# ###########################################################################
# ################ Start point2 2:30pm 20th Apr 2021  ########

# library(foreach)
# library(doParallel)

#Wins_all_location <- "Outputs/Wins_all_plasma.csv"
#Wins_all_plasma <- read_csv(Wins_all_location)
#Wins_all_plasma <- rename(Wins_all_plasma, ICC = IC, ICC_tot = IC_tot)

# Wap <- Wins_all
# #Wap <- mutate(Wap, S_tot = double(1))
# # CHUNK ND
# 
# ## Then regroup
# args <- commandArgs(TRUE)
# m1 <- as.numeric(args[1])
# m2 <- as.numeric(args[2])
# Wap <- Wap %>% filter(M >= m1 & M <= m2)
# 
# 
# 
# 
# #Wins_all_plasma$Q <- Wins_all_plasma$Ee + Wins_all_plasma$Ebind
# 
# 
# # S column means S_p
# # n_index <- read_csv("Dependencies/n_to_ICshell_Conversion.csv")
# #test <- filter(Wins_all_plasma_plot, Ef < 14.43)
# #### SO currently outputting the total NEEC rate and total resonance strength
# j<-1
# #registerDoParallel(15) 
# #foreach (i=1:length(Wap$ICC_tot)) %dopar% {
# for(i in 1:length(Wap$ICC_tot)) {
#   #if(  !is.na( Wap$B[i]) && !is.na( Wap$Ji_double[i]) && !is.na( Wap$Jf_double[i]) && !is.na( Wap$ICC[i]) &&  Wap$Type[i] != "E0" &&  Wap$Type[i] != "Large_L") {
#   #print(i)
#   #i <- 1
#   ##Replace
#   Wap$S_tot[i] <- compute_S(Wap$Jf_double[i], 1.5, Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i] )
#   #Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
#   Wap$N_neec_pw_Stot[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "N_neec_pw", "tot")
#   Wap$N_neec_pw[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "N_neec_pw", "p")
#   #}
#  
#     #Wins_all_plasma$shell[i] <- n_index$Shell[Wins_all_plasma$Occ[i]]
#   
#   ## Making S_n ... more accurate
# }
# 
# # plan(multisession, workers = 20)
# # #### DO it with pmap
# #    #Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , 1.5, Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i] ) 
# #    Wap$S_tot <- future_pmap(list(Wap$Jf_double, 1.5, Wap$Ji_double, Wap$ICC_tot, Wap$Ar, Wap$Q, Wap$Ebind ), compute_S)
# # #   #Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
# #   # Wap$N_neec_pw_Stot[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "N_neec_pw", "tot")
# #    #Wap$N_neec_pw[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "N_neec_pw", "p")
# Wins_all_location_save <- sprintf("Outputs/Chunk_Ni/Wins_all_plasma_%003.f_%003.f.csv", m1, m2)
# 
# #Wap$ICC_tot <- unlist(Wap$ICC_tot)
# #unlink(Wins_all_location)
# #if(!file.exists(Wins_all_location_save)){
# write.table(Wap, file=Wins_all_location_save, append=T, row.names=F, col.names=T,  sep=",")
# #}
# 
# 


# 
# j <- 1
# ### Rates Calculation for loop
# for(i in 1:length(Wap$S_tot)){ 
#   #i<-1
#   
#   print(i)
#   ## EBIT RATE
#   #AnN0$EBIT_Rate[i] <- compute_EBIT_RATE(AnN0$S[i], AnN0$Ee[i], V, nion, nelec, spread)
#   
#   j <- 1
#   #################### Needs J indexing properly
#   #for( j in 1:length(facility$Irradiance_Wcm2um2) ){
#   ## PLASMA RATE
#   if( is.na(Wap$S[i]) || is.na(Wap$Ee[i]) || is.na(Wap$Thf_numeric[i]) ) { 
#     Wap$N_neec_pw[i] = NA
#     Wap$Unity_Plasma_Rate[i] = NA
#     Wap$Ni[i] <- NA
#   }
#   else {
#     #Wins_all$N_neec_pw[i] <- compute_Plasma_Rate(Wins_all[i,], facility[j,])
#     Wap$N_neec_pw[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "N_neec_pw")
#     Wap$Unity_Plasma_Rate[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Rate")
#     Wap$Ni[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], "Nion")
#     
#   }
#   #}
#   
#   ##### compute rate for each facility
#   #Wins_all
#   
#   ## Charge state distribution rate
# }
# 
# 
# ### ADD SHELL TYPE FOR PLOTTING
# n_shell <- string_extract 
# 

### Save in Master_Rates file
  
