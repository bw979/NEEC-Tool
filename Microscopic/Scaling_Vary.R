library(rsconnect)
library(shiny)
library(tidyverse)
library(plotly)
library(RColorBrewer)


## set working directory as Microscopic/Scaling directory
#setwd("Microscopic/Scaling")

source("Wins.R")

## Read required data tables ... these are the newest tables for conversion of quantum labels
n_index <- read_csv("Microscopic/n_to_ICshell_Conversion.csv")
Ebinds <- read.csv("Dependencies/Ebinds/NIST_Ebinds.csv", header=TRUE, fill=TRUE)
#LS_Jat_Conv <- read.csv("Microscopic/elec_Jconv.csv")
## Add a dirac angular momentum quantum number variable

## This is actually read in Wins.R
#ICC_multipliers <- read.csv("Microscopic/ICC_Multipliers.csv")

## Read the Abinit data 
#### Contains all the wrangled Palffy data at this point
Abinit <- read.csv("Microscopic/Abinit_Apr2022.csv")

#Abinit <- filter(Abinit, GS==TRUE)

## Switch statement for divisor 
ICC_element <- ICC_Multipliers[which(ICC_Multipliers$Z==42),]
#ICC_ss <- ICC_element[which(ICC_element$subshell=="M5"),]

#Abinit  <- filter(Abinit, GS==TRUE)
Abinit <- mutate(Abinit, nh = double(1), nmax=double(1), Vi_q0=double(1), multiplier=double(1), a=double(1), b=double(1), Ratio=double(1))

for(i in 1:length(Abinit$q)){
  #print(i)
  Abinit$nh[i] <- ICC_element$nh[which(ICC_element$Occ == Abinit$Occ_after[i])]
  Abinit$nmax[i] <- ICC_element$nmax[which(ICC_element$Occ == Abinit$Occ_after[i])]
  Abinit$Vi_q0[i] <- ICC_element$Vi_q0[which(ICC_element$Occ == Abinit$Occ_after[i])]
}

########### CAUTION, MUST REMOVE THIS LINE AFTER PROPERLY SCALING A SHELL ##############
#Abinit <- filter(Abinit, n == 3)
########################################################################################
## ^^^^^^^^^^^^^^^ CAUTION ^^^^^^^^^^^^^^^^^^^^^^^^^^^ ######

multiplier <- function(a, b){
  m <- ((Abinit$Vi_keV[k]/Abinit$Vi_q0[k])^a)*((Abinit$nh[k]/Abinit$nmax[k])^b)
  S <- compute_S(8.5, Abinit$J_final[k], 10.5, Abinit$ICC[k]*m, Abinit$Ar[k], 4.85, Abinit$Vi_keV[k])
  val <- (S/Abinit$S_neec_beV[k]) 
  #rate_An <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$S_neec_beV[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i]
  
  return(val)
}



for(k in 1:length(Abinit$q)){
  
Abinit$ICC[k] <- Compute_ICC_singles("E2", Abinit$subshell[k], 42, 4.85, Abinit$Occ_after[k], "p", FALSE, 1, "")
#Abinit <- mutate(Abinit, nh = double(1), nmax=double(1), Vi_q0=double(1), multiplier=double(1), a=double(1), b=double(1), Ratio=double(1))
### Assign the multiplier variables
#Abinit <- Abinit
# multiplier <- function(a, b){
#   m <- ((Abinit$Vi_keV[k]/Abinit$Vi_q0[k])^a)*((Abinit$nh[k]/Abinit$nmax[k])^b)
#   S <- compute_S(8.5, Abinit$J_final[k], 10.5, Abinit$ICC[k]*m, Abinit$Ar[k], 4.85, Abinit$Vi_keV[k])
#   val <- (S/Abinit$S_neec_beV[k]) 
#   #rate_An <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$S_neec_beV[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i]
# 
#   return(val)
# }

###### Vary the power of the scaling fractions ATM JUST SETTING a=b=1
# a_vec <- c(seq(-2,4,1))
# b_vec <- c(seq(-2,4,1))
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

## Settings for all shells
 Abinit$a[k] <- 1
 Abinit$b[k] <- 1

m <- ((Abinit$Vi_keV[k]/Abinit$Vi_q0[k])^Abinit$a[k] )*((Abinit$nh[k]/Abinit$nmax[k])^Abinit$b[k] )
Abinit$multiplier[k] <- m
Abinit$Spdb[k] <- compute_S(8.5, Abinit$J_final[k], 10.5, Abinit$ICC[k], Abinit$Ar[k], 4.85, Abinit$Vi_keV[k]) * m
Abinit$Ratio[k] <- Abinit$Spdb[k] / Abinit$S_neec_beV[k]

}

Abinit <- mutate(Abinit, Eres_keV = Eres_eV/1000)

#dat
## Add atomic colours
cpal <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(13))

Abinit_GS <- filter(Abinit, GS==TRUE)
Residuals <- plot_ly(Abinit_GS, x =~Eres_keV, y =~q, z=~Ratio, type="scatter3d", mode="markers",
                     name = ~subshell, color = ~subshell, colors = cpal, size=1 ) %>%
  #add_segments(x = Abinit$Eres_keV[1], xend = Abinit$Eres_keV[1], y = Abinit$q[1], yend = Abinit$q[1], z=0, zend=Abinit$S_neec_beV[1], line=list(color = 'Black')) %>%
  #add_markers
  #add_paths(Abinit, x = ~Eres_keV, y = ~q, z = ~Ratio) %>%
  layout(
    scene=list(
      xaxis = list(title = "Eres (keV)" ),
      yaxis = list(title = "q"),
      zaxis = list(title = "Spdb / Sth",  exponentformat = 'E', dtick=10)
    )
  ) %>%
  config(mathjax = "cdn") #%>%
  #add_segments(x = Abinit_GS$Eres_keV[1], xend =Abinit_GS$Eres_keV[1], y = Abinit_GS$q[1], yend = Abinit_GS$q[1],z=0, zend=Abinit_GS$Ratio[1])
  
Residuals





############## START PLASMA RATE CALCULATION  ################
################## THIS HAS BEEN MOVED TO Macroscopic/FULL_Plasma_Macroscopic.R
# ne_set <- 1E24
# #Abinit <- mutate
# 
# ### RUN over the NEEC spectrum energies and decide on best Te for Ee = Te ####
# ### table for temps and rates
# Temps <- tibble(Te=(FLYCHK_Te*1E-3), Rtot=double(1), Rtot_psum = double(1), Rtot_PDBsum=double(1))
# for(j in 1:length(Temps$Te)){
# #### At this point S_sum is an order of magnitude too large
# #### Using Stot from alpha_tot = 4.63E5
# ######### Get average charge state #############
# ne <- ne_set ### e cm^-3
# #Te <- 1.441
# #<-22
# Te <- Temps$Te[j] ###keV
# Te_eV <- Te * 1E3
# Z_name <- sprintf("Dependencies/CSD/CSD_Ave_Data/Z%003.f_AVE_CS.csv", 42)
# Z_dat <- read_csv(Z_name)
# ##### Read the correct CSD fq file ####
# ## Round to nearest ne and Te
# #### test choose nearest val from Te and ne veactors
# # FLYCHK_ne <- c(1E12, 1E13, 1E14, 1E15, 1E16, 1E17, 1E18, 1E19, 1E20, 1E21, 1E22, 1E23, 1E24)
# # FLYCHK_Te <- c(0.5, 1, 1.5, 2, 5, 7, 10, 15, 23, 32, 52, 74, 100, 165, 235, 310, 390, 475, 655, 845, 1000, 1441, 1925, 2454, 3030, 3655, 4331, 5060, 5844, 6685, 7585, 8546, 10000, 20000, 50000, 100000)
# Te_ind <- which.min(abs(FLYCHK_Te - Te_eV))[1]
# Te_choice <- FLYCHK_Te[Te_ind]
# ne_val <- ne
# ne_ind <- which.min(abs(FLYCHK_ne - ne_val))[1]
# ne_choice <- FLYCHK_ne[ne_ind]
# CSD_File <- sprintf("Dependencies/CSD/CSD_fq_Data/CSD_%d.csv", 42)
# CSD <- read_csv(CSD_File)
# CSD_neTe <- filter(CSD, ne_cm3 == ne_choice, Te_eV == Te_choice)
# 
# #### We read the charge state of the CS chosen by the closest Te and ne selected by the laser params
# Z_ave <- filter(Z_dat, Te == Te_choice, ne == ne_choice)$CS
# #fq <- CSD_neTe$frac[which(CSD_neTe$q==(Wins_all$CS+1))]
# #if(Wins_all$Z >)
# 
# 
# ####### USING Spdb and S_palffy #######
# ###### Compare with plasma yields #####
# ## Literature Yield .... Te = 2 keV, ne = 1E23, Rneec = 5 s^-1
# ## S_tot Yield
# # ne <- 1E23 ### e cm^-3
# # Te <- 2 ###keV
# alpha_tot <- 4.63E5
# #compute_S <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){
# 
# ## Volume and number of ions
# Rp<- 40E-4 #cm
# Vp <- (4/3) * pi * (Rp^3) #cm^3
# ni <- (ne / Z_ave)
# f_iso <- 1E-5
# Ni <- ni * Vp * f_iso
# 
# 
# ### assign fq based on current {Te,ne} ###
# for(i in 1:length(Abinit$q)){
#   Abinit$fq[i] <- CSD_neTe$frac[which(CSD_neTe$q==(Abinit$q[i]))]
# }
# # for(i in 1:length(Mo$CS)){
# #   Mo$fq[i] <- CSD_neTe$frac[which(CSD_neTe$q==(Mo$CS[i]+1))]
# # }
# 
# 
# ### Calculate neec rate for each capture channel
# for(i in 1:length(Abinit$q)){
#   Abinit$R_neec_p[i] <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$S_neec_beV[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i]
#   Abinit$R_neec_pdb[i] <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$Spdb[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i]
#   
# }
# # for(i in 1:length(Mo$S)){
# #  Mo$Rneec_Pl[i] <- ne * MB(Mo$Ee[i]*1E3, 1E3*Te) *  sqrt(1 - (1/((1+(Mo$Ee[i]*1E3) / me))^2))*100*c   *  Mo$S[i] * 1E-24 * (2/(pi)) *  Mo$fq[i]
# # }
# R_neec_pTot <- sum(Abinit$R_neec_p)
# R_neec_PDBTot <-sum(Abinit$R_neec_pdb)
# 
# 
# ############# USING ALPHA TOT #############
# #### What is the best energy that gives the best sum rate
# # set parameters for sum rate
# #ne <- 1E24 ### e cm^-3
# #Te <- 1.441 ###keV
# # Alpha tot
# #Vi_chose <- mean(Abinit$Vi_keV)
# Vi_chose <-  mean(filter(ICC_Multipliers, Z==42)$Ebind)
# Q<-4.85
# alpha_tot <- 4.63E5
# Ee_chose <- Q-Vi_chose
# #E_max <- Q - mean(filter(ICC_Multipliers, Z==42)$Ebind)
# #E_max <- Q - max(Abinit$Vi_keV)
# S_beV <- compute_S(8.5, mean(Abinit$J_final), 10.5, alpha_tot, Abinit$Ar[10], Q, Vi_chose)
# Rp<- 40E-4 #cm
# 
# S <- 1E-24 * S_beV #cm^2 eV
# F_ETe <- MB(1E3*Ee_chose, 1E3*Te)
# v_max <- sqrt(1 - (1/((1+(Ee_chose*1E3) / me))^2))  #c
# v_max <- 100*c*v_max #cm s^-1
# 
# ## Volume and number of ions
# Vp <- (4/3) * pi * (Rp^3) #cm^3
# ni <- (ne / Z_ave)
# f_iso <- 1E-5
# Ni <- ni * Vp * f_iso
# 
# ## Flux and cross section and adjustables
# flux <- ne *  v_max * F_ETe
# N_res <- 1
# sigma_res <- S * (2/(pi))
# 
# R_neec_tot <- S * (2/(pi)) * F_ETe * N_res * ne * v_max
# 
# Temps$Rtot[j] <- R_neec_tot
# Temps$Rtot_psum[j] <- R_neec_pTot
# Temps$Rtot_PDBsum[j] <- R_neec_PDBTot
# 
# # R_neec_tot
# # R_neec_pTot
# 
# }


# Temps$R_Alpha_tot <- Temps$Rtot
# Temps$Rtot_PalffyWu[j] <- Temps$Rtot_psum
# Temps$Rtot_Wallis <- Temps$Rtot_PDBsum

#### PDB S Yield



##### Look at average charge states over Te and ne   #################


#### Configuration reader code 



# ### PLOTTING rates ########
# pE20 <- plot_ly(Temps, x = ~Te, y = ~Rtot, type = 'scatter', mode='markers+lines', showlegend=TRUE, name="R_Alpha_tot", color="red") %>% #, marker = list(color = label_colours)) %>%
#   layout(yaxis = list(  title = TeX("R_{neec} (s^{-1})"), exponentformat="E"), #range = c(1E-10, max(wap_filtered$S))),
#          xaxis = list(tickangle = -45, title="Temp [keV]", range=c(0.1,30)) #%>%
#          #legend=list(name=~shell, showlegend=TRUE) %>%
# 
#          # yaxis2 = list(
#          #   tickfont = list(color = "red"),
#          #   overlaying = "y",
#          #   side = "right",
#          #   title = "MB" ,
#          #   showlegend = F)
# 
#   )  %>%
#   add_trace(x= ~Te, y= ~Rtot_psum, type='scatter', mode="markers+lines", name="Rtot_PalffyWu", color="green") %>%
#   add_trace(x= ~Te, y= ~Rtot_PDBsum, type='scatter', mode="markers+lines", name = "Rtot_Wallis", color="yellow") %>%
#   config(mathjax = "cdn")
# 
# pE20
# 
# fig <- subplot(pE19, pE20, pE21, pE22, pE23, pE24, nrows = 6)
# fig

  #add_lines(x = ~Evals, y = ~MBvals, type = 'scatter', mode = 'lines', color="#DE3163", showlegend=T, name="F(E)")


###### PLOT S_pdb / S_th
#Abinit <- Abinit_GS

# Abinit <- mutate(Abinit, Eres_keV = Eres_eV/1000)
# 
# Residuals <- plot_ly(Abinit, x =~Eres_keV, y =~q, z=~Ratio, type="scatter3d", mode="markers",
#                      name = ~subshell, color = ~subshell, marker = list(size=2.51) ) %>%
#   #add_segments(x = Abinit$Eres_keV[1], xend = Abinit$Eres_keV[1], y = Abinit$q[1], yend = Abinit$q[1], z=0, zend=Abinit$S_neec_beV[1], line=list(color = 'Black')) %>%
#   #add_markers
#   #add_paths(Abinit, x = ~Eres_keV, y = ~q, z = ~Ratio) %>%
#   layout(
#     scene=list(
#       xaxis = list(title = "Eres_keV" ),
#       yaxis = list(title = "q"),
#       zaxis = list(title = "Spdb / Sth", type='log', exponentformat = 'E', range=c(0,2))
#     )
#   )
# Residuals

wap <- read_csv("Outputs/MASTER/All/All_RatesCalcd_Viva.csv")
#wap <- bind_rows(wap, additions)

#### CHecking out an X-NEEC Method for 73Ge
data <- filter(wap, AX=="73GE", Q==13.2845)
Q <- 13.2845
#Vnew <- filter(ICC_Multipliers, Z==data$Z[1])$Vi_q0[1]
Vnew <-  Q - (1/100)*Q

#Qnew <- Vnew + (1/10)*Vnew
#Qnew <- Q - (1/100)*Q


#Compute_ICC_singles("E2", Abinit$subshell[k], 42, 4.85, Abinit$Occ_after[k], "p", FALSE, 1, "")
Compute_ICC_singles("E2", "K" , data$Z[1], Qnew, 1, "p", FALSE, NA, NA)
S_Xneec <- compute_S(data$Jf_double[1], data$Jat[1], data$Ji_double[1], 299, data$Ar[1], data$Q[1], Vnew)









#Abinit$Spdb <- compute_S(8.5, Abinit$J_final[i], 10.5, Abinit$ICC[i], Abinit$Ar[i], 4.85, Abinit$Vi_keV[i]) 



