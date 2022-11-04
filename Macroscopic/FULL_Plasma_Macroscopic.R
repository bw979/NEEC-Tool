library(rsconnect)
library(shiny)
library(tidyverse)
library(plotly)

#### THIS ONE COMPARES THE MACROSCOPIC RATE TO THE PALFFY PUBLISHED RATES LOOKING AT ICC SCALING FROM DIFFERENT PERSPECTIVES

## set working directory as Microscopic/Scaling directory
#setwd("Microscopic/Scaling")

source("Wins.R")

## Read required data tables ... these are the newest tables for conversion of quantum labels
n_index <- read_csv("Microscopic/n_to_ICshell_Conversion.csv")
Ebinds <- read.csv("Dependencies/Ebinds/NIST_Ebinds.csv", header=TRUE, fill=TRUE)
Abinit <- read.csv("Microscopic/Abinit_Apr2022.csv")
additions <- read_csv("Mo_Additions_Rates2.csv")
wap <- read_csv("Outputs/MASTER/All/All_RatesCalcd_Viva.csv")
wap <- bind_rows(wap, additions)
ICC_element <- ICC_Multipliers[which(ICC_Multipliers$Z==42),]




### Need to add in the wap data to the abinit data.frame
Abinit <- mutate(Abinit, S_wap = double(1), R_neec_wap = double(1))
wap_f2 <- filter(wap, AX == "93MO", Q < 4.86) %>%
  filter(!is.na(S))
wap_f2 <- mutate(wap_f2, S2 = double(1), ICC2 = double(1))
## Check S calculation
#wap_f2$ICC2 <-  pmap(list(wap_f2$Type, wap_f2$Z, wap_f2$Q, "p"), Compute_ICC_singles )

for(i in 1:length(wap_f2$S)){
  #Type, ss,  Zin, Q, n, OUTPUT, MIX, MR, MIX_type)
  ### Recalculate unscaled S (in wap it has been scaled by the gamma detection factor)
  wap_f2$ICC2[i] <- Compute_ICC_singles(wap_f2$Type[i], wap_f2$subshell[i], wap_f2$Z[i], wap_f2$Q[i], wap_f2$Occ[i], "p", FALSE, NA, NA)
  wap_f2$S2[i] <- compute_S(8.5, wap_f2$Jat[i], 10.5, wap_f2$ICC2[i], wap_f2$Ar[10], wap_f2$Q[i], wap_f2$Ebind[i])
}

#for viewing the data
wap_f2_view <- select(wap_f2, -MIX_type, -B2, -B2type, -FL)

for(i in 1:length(Abinit$q)){
 print(i)
if(Abinit$GS[i]) { 
 if(is.null(wap_f2$S2[which(wap_f2$Occ==Abinit$Occ_after[i])]  )) {Abinit$S_wap[i] <- NA}
 if(is.na(wap_f2$S2[which(wap_f2$Occ==Abinit$Occ_after[i])]  )) {Abinit$S_wap[i] <- NA}
 Abinit$S_wap[i] <- wap_f2$S2[which(wap_f2$Occ==Abinit$Occ_after[i])]
} else {
  Abinit$S_wap[i] <- NA
}
}

############# SORTING OUT THE MULTIPLIERS  ###################
## Look at a configuration averaging of the NEEC-X channels


Abinit <- mutate(Abinit, nh = double(1), nmax=double(1), Vi_q0=double(1), multiplier=double(1), a=double(1), b=double(1), Ratio=double(1))
## Assign scaling variables
for(i in 1:length(Abinit$q)){
  #print(i)
  Abinit$nh[i] <- ICC_element$nh[which(ICC_element$Occ == Abinit$Occ_after[i])]
  Abinit$nmax[i] <- ICC_element$nmax[which(ICC_element$Occ == Abinit$Occ_after[i])]
  Abinit$Vi_q0[i] <- ICC_element$Vi_q0[which(ICC_element$Occ == Abinit$Occ_after[i])]
}

for(k in 1:length(Abinit$q)){
  ## Assign an unscaled ICC
  Abinit$ICC[k] <- Compute_ICC_singles("E2", Abinit$subshell[k], 42, 4.85, Abinit$Occ_after[k], "p", FALSE, 1, "")
 
  multiplier <- function(a, b){
    m <- ((Abinit$Vi_keV[k]/Abinit$Vi_q0[k])^a)*((Abinit$nh[k]/Abinit$nmax[k])^b)
    S <- compute_S(8.5, Abinit$J_final[k], 10.5, Abinit$ICC[k]*m, Abinit$Ar[k], 4.85, Abinit$Vi_keV[k])
    val <- (S/Abinit$S_neec_beV[k]) 
    return(val)
  }
  
  # a_vec <- c(seq(-2,4,0.1))
  # b_vec <- c(seq(-2,4,0.1))
  # dl <- length(a_vec)*length(b_vec)
  # 
  # #multi <- tibble(a=double(dl), b=double(dl), Ratio=double(dl))
  # multi <- crossing(a_vec, b_vec) %>%
  #   mutate(Ratio=double(dl)) %>%
  #   rename(a=a_vec, b=b_vec)
  # 
  # 
  # for(i in 1:dl){
  #   multi$Ratio[i] <- multiplier(multi$a[i], multi$b[i])
  # }           
  # 
  # Abinit$a[k] <-multi[which.min(abs(multi$Ratio - 1))[1],]$a
  # Abinit$b[k] <- multi[which.min(abs(multi$Ratio - 1))[1],]$b
  
  #Abinit$a[k] <- Abinit$nh[k] - Abinit$nmax[k]
  #Abinit$b[k] <- Abinit$nmax[k] - Abinit$nh[k]
  
  Abinit$a[k] <- 1
  Abinit$b[k] <- 1
  
  
  ### Try scaling
  
  
  
  m <- ((Abinit$Vi_keV[k]/Abinit$Vi_q0[k])^Abinit$a[k] )*((Abinit$nh[k]/(Abinit$nmax[k]))^Abinit$b[k] )
  Abinit$multiplier[k] <- m
  #### CURRENTLY UNSCALED Spdb
  Abinit$Spdb[k] <- compute_S(8.5, Abinit$J_final[k], 10.5, Abinit$ICC[k], Abinit$Ar[k], 4.85, Abinit$Vi_keV[k])#*m
  Abinit$Ratio[k] <- Abinit$Spdb[k] / Abinit$S_neec_beV[k]
  
}



############## START PLASMA RATE CALCULATION  ################
ne_set <- 1E20

### RUN over the NEEC spectrum energies and decide on best Te for Ee = Te ####
### table for temps and rates
Temps <- tibble(Te=(FLYCHK_Te*1E-3), Rtot=double(1), Rtot_psum = double(1), Rtot_PDBsum=double(1), Rtot_wap=double(1))
for(j in 1:length(Temps$Te)){
  #### At this point S_sum is an order of magnitude too large
  #### Using Stot from alpha_tot = 4.63E5
  ######### Get average charge state #############
  ne <- ne_set ### e cm^-3
  #Te <- 1.441
  #<-22
  Te <- Temps$Te[j] ###keV
  Te_eV <- Te * 1E3
  Z_name <- sprintf("Dependencies/CSD/CSD_Ave_Data/Z%003.f_AVE_CS.csv", 42)
  Z_dat <- read_csv(Z_name)
  ##### Read the correct CSD fq file ####
  ## Round to nearest ne and Te
  #### test choose nearest val from Te and ne veactors
  # FLYCHK_ne <- c(1E12, 1E13, 1E14, 1E15, 1E16, 1E17, 1E18, 1E19, 1E20, 1E21, 1E22, 1E23, 1E24)
  # FLYCHK_Te <- c(0.5, 1, 1.5, 2, 5, 7, 10, 15, 23, 32, 52, 74, 100, 165, 235, 310, 390, 475, 655, 845, 1000, 1441, 1925, 2454, 3030, 3655, 4331, 5060, 5844, 6685, 7585, 8546, 10000, 20000, 50000, 100000)
  Te_ind <- which.min(abs(FLYCHK_Te - Te_eV))[1]
  Te_choice <- FLYCHK_Te[Te_ind]
  ne_val <- ne
  ne_ind <- which.min(abs(FLYCHK_ne - ne_val))[1]
  ne_choice <- FLYCHK_ne[ne_ind]
  CSD_File <- sprintf("Dependencies/CSD/CSD_fq_Data/CSD_%d.csv", 42)
  CSD <- read_csv(CSD_File)
  CSD_neTe <- filter(CSD, ne_cm3 == ne_choice, Te_eV == Te_choice)
  
  #### We read the charge state of the CS chosen by the closest Te and ne selected by the laser params
  Z_ave <- filter(Z_dat, Te == Te_choice, ne == ne_choice)$CS
  #fq <- CSD_neTe$frac[which(CSD_neTe$q==(Wins_all$CS+1))]
  #if(Wins_all$Z >)
  
  
  ####### USING Spdb and S_palffy #######
  ###### Compare with plasma yields #####
  ## Literature Yield .... Te = 2 keV, ne = 1E23, Rneec = 5 s^-1
  ## S_tot Yield
  # ne <- 1E23 ### e cm^-3
  # Te <- 2 ###keV
  alpha_tot <- 4.63E5
  #compute_S <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){
  
  ## Volume and number of ions
  Rp <- 40E-4 #cm
  Vp <- (4/3) * pi * (Rp^3) #cm^3
  ni <- (ne / Z_ave)
  f_iso <- 1E-5
  Ni <- ni * Vp * f_iso
  
  
  ### assign fq based on current {Te,ne} ###
  for(i in 1:length(Abinit$q)){
    Abinit$fq[i] <- CSD_neTe$frac[which(CSD_neTe$q==(Abinit$q[i]))]
  }
  # for(i in 1:length(Mo$CS)){
  #   Mo$fq[i] <- CSD_neTe$frac[which(CSD_neTe$q==(Mo$CS[i]+1))]
  # }
  
  
  ### Calculate neec rate for each capture channel
  for(i in 1:length(Abinit$q)){
    Abinit$R_neec_p[i] <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$S_neec_beV[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i]
    Abinit$R_neec_pdb[i] <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$Spdb[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i]
    #Abinit$R_neec_wap[i] <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$S_wap[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i]
  }
  Abinit_GS <- filter(Abinit, GS==TRUE)
  for(i in 1:length(Abinit_GS$q)){
    Abinit_GS$R_neec_wap[i] <- ne * MB(Abinit_GS$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit_GS$Eres_eV[i]) / me))^2))*100*c   *  Abinit_GS$S_wap[i] * 1E-24 * (2/(pi)) *  Abinit_GS$fq[i] * Abinit_GS$multiplier[i]
  }  
  
  # for(i in 1:length(Mo$S)){
  #  Mo$Rneec_Pl[i] <- ne * MB(Mo$Ee[i]*1E3, 1E3*Te) *  sqrt(1 - (1/((1+(Mo$Ee[i]*1E3) / me))^2))*100*c   *  Mo$S[i] * 1E-24 * (2/(pi)) *  Mo$fq[i]
  # }
  R_NEEC_pTot <- sum(Abinit$R_neec_p)
  R_NEEC_PDBTot <-sum(Abinit$R_neec_pdb)
  R_NEEC_wapTot <- sum(Abinit_GS$R_neec_wap)
  
  
  ############# USING ALPHA TOT #############
  #### What is the best energy that gives the best sum rate
  # set parameters for sum rate
  #ne <- 1E24 ### e cm^-3
  #Te <- 1.441 ###keV
  # Alpha tot
  #Vi_chose <- mean(Abinit$Vi_keV)
  #Vi_chose <-  mean(filter(ICC_Multipliers, Z==42)$Ebind)
  Vi_chose <- mean(wap_f2$Ebind)
  Q<-4.85
  alpha_tot <- 4.63E5
  Ee_chose <- Q-Vi_chose
  #E_max <- Q - mean(filter(ICC_Multipliers, Z==42)$Ebind)
  #E_max <- Q - max(Abinit$Vi_keV)
  S_beV <- compute_S(8.5, mean(Abinit$J_final), 10.5, alpha_tot, Abinit$Ar[10], Q, Vi_chose)
  Rp<- 40E-4 #cm
  
  S <- 1E-24 * S_beV #cm^2 eV
  F_ETe <- MB(1E3*Ee_chose, 1E3*Te)
  v_max <- sqrt(1 - (1/((1+(Ee_chose*1E3) / me))^2))  #c
  v_max <- 100*c*v_max #cm s^-1
  
  ## Volume and number of ions
  Vp <- (4/3) * pi * (Rp^3) #cm^3
  ni <- (ne / Z_ave)
  f_iso <- 1E-5
  Ni <- ni * Vp * f_iso
  
  ## Flux and cross section and adjustables
  flux <- ne *  v_max * F_ETe
  N_res <- 1
  sigma_res <- S * (2/(pi))
  
  R_NEEC_tot <- S * (2/(pi)) * F_ETe * N_res * ne * v_max
  
  Temps$Rtot[j] <- R_NEEC_tot
  Temps$Rtot_psum[j] <- R_NEEC_pTot
  Temps$Rtot_PDBsum[j] <- R_NEEC_PDBTot
  Temps$Rtot_wap[j] <- R_NEEC_wapTot
  
  # R_neec_tot
  # R_neec_pTot
  
}

# yaxis2 = list(
#   tickfont = list(color = "red"),
#   overlaying = "y",
#   side = "right",
#   title = "MB" ,
#   showlegend = F)

  
  
  
## Graphs of rates for different methods
### PLOTTING ########
pE20 <- plot_ly(Temps, x = ~Te, y = ~Rtot, type = 'scatter', mode='markers+lines', showlegend=TRUE, name=TeX("\\alpha_{Tot}"), line=list(color="orange"), marker=list(color="orange")) %>% #, marker = list(color = label_colours)) %>%
  layout(yaxis = list(  title = TeX("R_{NEEC} (Ion^{-1}s^{-1})"), exponentformat="E"), #range = c(1E-10, max(wap_filtered$S))),
         xaxis = list(tickangle = -45, title = TeX("T_e (keV)"), range=c(0.1,30)) #%>%
         #legend=list(name=~shell, showlegend=TRUE) %>%
         
  )  %>%
  add_trace(x= ~Te, y= ~Rtot_psum, type='scatter', mode="markers+lines", name="[PalffyWu] Ab Initio Data", line=list(color="green"), marker=list(color="green")) %>%
  add_trace(x= ~Te, y= ~Rtot_PDBsum, type='scatter', mode="markers+lines", name = "AbInit Energies and CSF's PDB Unscaled", line=list(color="red"), marker=list(color="red")) %>%
  #add_trace(x= ~Te, y= ~Rtot_PDBsum, type='scatter', mode="markers+lines", name = "Rtot_pdata_allUnscaled", line=list(color="red"), marker=list(color="red")) %>%
  add_trace(x= ~Te, y= ~Rtot_wap, type='scatter', mode="markers+lines", name = "PDB - Scaling a=b=1", line=list(color="yellow", marker=list(color="yellow"))) %>%
  config(mathjax = "cdn") %>%
  layout(legend= list(x=0.6, y=0.9))

pE20


######### Plot full abinit data set
# Abinit$Eres_keV <- Abinit$Eres_eV * 1E-3
# 
# plot_ly(Abinit, x =~Eres_keV, y = ~S_neec_beV, type = 'bar', showlegend=TRUE, name=~shell)
# 
# 
# 





