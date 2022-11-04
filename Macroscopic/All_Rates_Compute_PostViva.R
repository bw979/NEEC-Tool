library(tidyverse)
library(ggplot2)
library(stats)
library(grid)
library(plotly)

source("Wins.R")

#### SEPT 2022 - Currently this evaluates the alpha_tot FOM for all EBIT and PLASMAS at ELINP
## Generates Wap data for all NEECable An channels


##########################
######## PLASMA ##########
##########################
FLYCHK_ne <- c(1E12, 1E13, 1E14, 1E15, 1E16, 1E17, 1E18, 1E19, 1E20, 1E21, 1E22, 1E23, 1E24)
FLYCHK_Te <- c(0.5, 1, 1.5, 2, 5, 7, 10, 15, 23, 32, 52, 74, 100, 165, 235, 310, 390, 475, 655, 845, 1000, 1441, 1925, 2454, 3030, 3655, 4331, 5060, 5844, 6685, 7585, 8546, 10000, 20000, 50000, 100000)
##For the subshell column
ICC_Multipliers <- read_csv("Dependencies/ICC_Multipliers.csv") 
#ICC_Multipliers <-  arrange(ICC_Multipliers, Z, Occ)
atomic_colours <- read_csv("Dependencies/Atomic_Colours.csv")


##############################################################
############# FULL ALPHA_tot CALC - 1 CALC PER Q #############
##############################################################
#### USE THIS BIT IF WE HAVE ALREADY HAD A RUN OF RATES
Wap_N0 <- read_csv("Outputs/MASTER/All/AnN0_all.csv")
Wap_Ni <- read_csv("Outputs/MASTER/All/AnNi_all.csv")
Wap <- bind_rows(Wap_N0, Wap_Ni)
additions <- read_csv("Mo_Additions_Rates2.csv")
#additions2 <- read_csv("Rb_Addition.csv")
Wap <- bind_rows(Wap, additions)




### MAY NEED SOME WRANGLING USING THIS ACTUALLY
#Wap <- read_csv("Outputs/MASTER/All/All_RatesCalcd.csv")
#Wap_PV <- read_csv("Outputs/MASTER/All/All_RatesCalcd_Viva.csv")
#Wap_PV <- bind_rows(Wap_PV, additions)
#rm(Wap_N0)
#rm(Wap_Ni)
########## Wrangle and save the full database at this point ande put it in final structure ###############
##### ENSURE YOU HAVE IRONED OUT SPECIAL CASE B Values in the GAMMA DATABASE
# ## Add mixing columns and ICC-tot
# Wap <- Wap %>%
#   mutate(Jf_double=double(1), Ji_double=double(1), MIX_type=character(1), ICC_tot=double(1), Stot=double(1),
#          Gamma_Thalf_eV=double(1), Gamma_eV=double(1), ICC_sum=double(1), Npw_Pl_S = double(1), S_sum=double(1), Npw_SEBIT=double(1) )
# 
# #### FINAL DATABASE STRUCTURE
# Wap <- select(Wap, AX, M, El, Z, N, Ei, Ef, Ji, Ji_double, Jf, Jf_double, Type, Thi, Thi_numeric, Thf, Thf_numeric,
#                Occ, CS, At_GS_Config, subshell, Jat, Ebind, Q, Egam, Efeed, FL, Btype, B, MR, MIX_type, B2type, B2, RI, Ar, Gamma_Thalf_eV, Gamma_eV,
#               ICC, ICC_tot, ICC_sum, S, Stot, S_sum, Ee, Ebeam, -Rate, Npw_Pl_S, Npw_SEBIT)

#does it find 83Mo - YES
#wap_f2 <- filter(Wap, AX == "93MO", Q < 4.86)
#wap_f3 <- filter(Wap, AX == "84RB", Q < 4)


### Arrange so can make FOM DATABASE
Wap <- arrange(Wap, M, Z, Occ )

#Wap <- select(Wap, -Rneec_Pl, Rneec_SEBIT)
###
Wap$Q <- Wap$Ef - Wap$Ei
Wap$Ee <- Wap$Q - Wap$Ebind
Wap <- filter(Wap, Ee>0, Z>11, Z<104) %>%
  mutate(B_WE = FALSE)

### Replace NA Bvalues with Bvalue = 1 with lapply

B_NA_Identify <- function(value){
  if(is.na(value)){
    return(TRUE)
  }
  else{ return(FALSE)}
}

B_NAto1 <- function(value){
  if(is.na(value)){
    return(1)
  }
  else {return(value)}
}

#### NEED to bind together Wap_PV and Wap_N0 and Ni first
Wap$B_WE <- lapply(Wap$B, B_NA_Identify)
Wap$B <- lapply(Wap$B, B_NAto1)
Wap$B <- as.double(Wap$B)

# #### FINAL DATABASE STRUCTURE from AnN0_all + AnNi_all
Wap <- mutate(Wap, Ji_double=double(1), Jf_double=double(1), MIX_type = character(1), Gamma_Thalf_eV=double(1), Gamma_eV=double(1), 
              ICC_tot = double(1), ICC_sum = double(1), m=double(1), ICC_Vac=double(1), ICC_p=double(1), S_Vac=double(1), Type_Good = FALSE, 
              Stot=double(1), Rate_Pl_tot=double(1), Rate_EBIT_tot=double(1), Ebind_mean=double(1), Jat_mean=double(1), Sp=double(1)) 

Wap <- select(Wap, AX, M, El, Z, N, Ei, Ef, Ji, Ji_double, Ji_double, Jf, Jf_double, Type, Type_Good, Thi, Thi_numeric, Thf, Thf_numeric,
               Occ, CS, At_GS_Config, subshell, Jat, Jat_mean, Ebind, Ebind_mean, Q, Ee, Ebeam, Egam, Efeed, FL, Btype, B, MR, MIX_type, B2type, B2, B_WE, RI, Ar, Gamma_Thalf_eV, Gamma_eV,
              ICC_p, ICC_tot, -ICC_sum, Sp, Stot, m, Ee, -Rate, Rate_Pl_tot, Rate_EBIT_tot) %>%
  filter(Type != "E0")


## Useful database of ALL the resonance channels
Wap_ALL <- Wap

### STILL NOT SAVED IT YET NEEDS RUNNING THROUGH TO CALC ICCp's
Wap_ALL_location_save <- "Wap_ALL_BWE_CorrectionsSubmission.csv"
## Save Wap_ALL database
if(!file.exists(Wap_ALL_location_save)){
  write.table(Wap_ALL, file=Wap_ALL_location_save, append=T, row.names=F, col.names=T,  sep=",")
}

Wap <- distinct(Wap, AX, Q, .keep_all = T)

Wap$Btype <- str_remove_all(Wap$Btype, "[BW\\=]")
Wap$B2type <- str_remove_all(Wap$Btype, "[BW\\=]")


####### NEED
##### PRESENT... WANT TO calculate ICC's and do some initial rate calcs
for(i in 1:nrow(Wap) ){
#for(i in 1:300 ){
#i<-80
  # # ### some extra wrangling to add subshell column (xray notation)
  # # #Wap$subshell[i] <- filter(ICC_Multipliers, Z==Wap$Z[i], Occ==Wap$Occ[i])$subshell
  # # Wap$shell[i] <- str_remove_all(Wap$subshell[i], "\\d")
  # # Wap$shell_colour[i] <- atomic_colours$colour2[which(atomic_colours$subshell==Wap$subshell[i])[1]] 
  # # 
  
  #### Fully Ionised Temperature Parameters
  #Z_name <- sprintf("Dependencies/CSD/CSD_Ave_Data/Z%003.f_AVE_CS.csv", Wap$Z[i])
  #Z_dat <- read_csv(Z_name)
  ## Find the min temp at that density
 # Te_min <- filter(Z_dat, ne==1E24, CS > (Wap$Z[i]-0.5))$Te[1]
  
  ## assign Te_min
  #if(!is.na(Te_min)){
  #  Wap$Te_Min[i] <- Te_min 
  # Wap$Etot[i] <- sum(filter(ICC_Multipliers, Z==Wap$Z[i])$Ebind)
  #} else {
    Wap$Etot[i] <- sum(filter(ICC_Multipliers, Z==Wap$Z[i])$Ebind)
    Wap$Te_min[i] <- Wap$Etot[i]
      
  #}
  
  #### Make sure the correct B is assigned
  if(is.na(Wap$Btype[i])){
    Wap$Type_Good[i] <- NA
  } else if(((Wap$Type[i] == Wap$Btype[i]) && !is.na(Wap$Btype[i])) || Wap$B_WE[i] == TRUE){
    Wap$Type_Good[i] = TRUE
  }
  
  ### Assign the B2 if necessary
  if(Wap$Type_Good[i] == FALSE &&  !is.na(Wap$Btype[i]) )
  {
    if(Wap$Type[i] == Wap$B2type[i]) {
      Wap$B[i] <- Wap$B2[i]
    } else {
      Wap$B_WE[i] <- TRUE
      Wap$B[i] <- 1
    }  
      
      
      
    # } else if(Wap$Type[i] == "M1") {
    #   Wap$B[i] <- Wap$B2[i]
    # } else {
    #   Wap$B_WE[i] <- TRUE
    #   Wap$B[i] <- 1
    # }
  }
  

  ### Create the mixing type column eg. "M1+E2"
  if( isFALSE(Wap$B_WE[i]) && Wap$Btype[i] =="BM1W" && Wap$MR[i]!=0 && !is.na(Wap$MR[i])  ){
#  if( !is.na(Wap$Btype[i]) && !is.na(Wap$B2type[i]) && Wap$Btype[i] =="BM1W" && Wap$B2type[i]=="BE2W=" && Wap$MR[i]!=0 && !is.na(Wap$MR[i])  ){
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
  #i<-1
  #### Evaluate Mean ionisation energy from ICC_Multipliers
  ICC_element <- filter(ICC_Multipliers, Z == Wap$Z[i])
  Wap$Ebind_mean[i] <- mean(filter(ICC_element, Occ>=Wap$Occ[i])$Ebind)
  Wap$Jat_mean[i] <- mean(filter(ICC_element, Occ>=Wap$Occ[i])$Jat)
  Wap$subshell[i] <- ICC_element$subshell[which(ICC_element$Occ == Wap$Occ[i])]
  
  ##### CALCULATE Ar, ICC and Stot
  #i <- 5000
  if(   !is.na(Wap$Ji_double[i]) && !is.na(Wap$Jf_double[i]) && !is.na(Wap$Type[i]) && Wap$Type[i] != "E0" && Wap$Type[i] != "Large_L") {
    # ######Compute_Ar... reuturns in SI [s^-1] ... also trying ICC_tot
    # #if its mixed
    # if(  (Wap$MR[i]!=0) && (!is.na(Wap$MR[i])) ){
    #   Wap$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i], TRUE, Wap$MR[i] )
    #   ## If its "M1 + E2" mixed
    #   if( !is.na(Wap$MIX_type[i]) && Wap$MIX_type[i] == "M1+E2"  ){
    #     try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot", TRUE,  Wap$MR[i], Wap$MIX_type[i])  )
    #     try(Wap$ICCp[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", TRUE, Wap$MR[i], Wap$MIX_type[i])  )
    #     # try(Wap$ICC_sum[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "SUMp", TRUE, Wap$MR[i], Wap$MIX_type[i])  )
    #   }
    # } else {
    
      Wap$m[i] <- ICC_element$multiplier[which(ICC_element$Occ == Wap$Occ[i])]
    
    
      Wap$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i], FALSE, Wap$MR[i] )
      # try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot") )
      # try(Wap$ICC_p[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p"))
      try(Wap$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$subshell[i], Wap$Z[i], Wap$Q[i], 1, "tot", FALSE, 1, "") )
      try(Wap$ICC_p[i] <- Compute_ICC_singles(Wap$Type[i], Wap$subshell[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", FALSE, 1, "") * Wap$m[i])
      try(Wap$ICC_frac[i] <- Compute_ICC_singles(Wap$Type[i], Wap$subshell[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", FALSE, 1, ""))
      ## The available ICC for capture ... must combine this with an effective energy
      Wap$ICC_Vac[i] <- Wap$ICC_tot[i] - Wap$ICC_frac[i] + Wap$ICC_p[i]
      # try(Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot"))
      # try(Wap$ICC_tot[i] <- Compute_ICC(Wap[i,], "tot"))
      #try(Wap$ICC_sum[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "SUMp", FALSE, 0, "")  )
    # }
  }
  #print(i)
  
  
  #   # #Compute_S Jd, Jat,Ji,ICC, Ar, Etrans_keV, Ebind_keV
  #   # Wap$S[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
  #   # Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
  # 
  # 
  # ### Total depletion level width
  Wap$Gamma_Thalf_eV[i] <- (hbar/e) * (log(2) / Wap$Thf_numeric[i])
  Wap$Gamma_eV[i] <-  (hbar/e) * (Wap$Ar[i]*(1+Wap$ICC_p[i]))
  #Wap$Gamma_Scale_Factor[i] <- 1/(1+Wap$ICC_sum[i])

  #### RATES ####

  #j <- which(facility$Facility == "ELI-NP")   #...ELI NP
  
    ### COMPUTE THE RESONANCE STRENGTH
    Wap$S_Vac[i] <-   compute_S(Wap$Jf_double[i], Wap$Jat_mean[i], Wap$Ji_double[i], Wap$ICC_Vac[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind_mean[i] )
    Wap$Stot[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat_mean[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind_mean[i] )
    Wap$Sp[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_p[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i]) 
    # Wap$S_sum[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_sum[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
    print(i)
    
  ne <- 1E24 #cm^-3
  Stot <- 1E-24 * Wap$Stot[i] #cm^2 eV
  #Sp <- 1E-24 * Wap$Sp[i] #cm^2 eV
  ## Choose the  Temp at the mean impact energy
  Ee_chose <- Wap$Q[i] - Wap$Ebind_mean[i]
  Te <- Ee_chose
  F_ETe <- MB(1E3*Ee_chose, 1E3*Te)
  v_chose <- sqrt(1 - (1/((1+(Ee_chose*1E3) / me))^2))  #c
  v_chose <- 100*c*v_chose #cm s^-1
  #try(Wap$Rate_Pl_tot[i] <- compute_Plasma_Rate(Wap[i,], facility[j,], 100E-12, "Te", (1000*facility$Te_Beg_keV[j]), "Rate", "p"))
  
  
  ##### BOTH OUTPUTTING PER ION CURRENTLY
  try( Wap$Rate_Pl_tot[i] <- ne * v_chose  * F_ETe  * Stot * (2/(pi)) )   ### Ion^-1 s^-1
  try(Wap$Rate_EBIT_tot[i] <- compute_EBIT_RATE(Wap$Sp[i], Wap$Ee[i], Wap$Gamma_eV[i], "Rate") )  ##s^-1
  # try(Wap$Npw_Pl_S[i] <- compute_Plasma_Rate(Wap[i,], facility[j,] , 100E-12, "Te", (1000*facility$Te_Beg_keV[j]), "N_neec_pw", "p"))
  # try(Wap$Npw_SEBIT[i] <- compute_EBIT_RATE(Wap$S[i], Wap$Ee[i], Wap$Gamma_eV[i], "N_neec_pw") )
  # 
}

# #Wins_all_location_save <- sprintf("Outputs/MASTER/Chunk/Wins_all_plasma_%003.f_%003.f.csv", m1, m2)

#### SUPPOSED TO SAY NO MIXING
Wins_all_location_save <- "Wap_AllQ3_PlasmaEBIT_NoMixing_CorrectedSubmission2022.csv"
###### DONT WRITE FOR NOW THANKS ###
###################################
Wap$B_WE <- as.character(Wap$B_WE)

if(!file.exists(Wins_all_location_save)){
  write.table(Wap, file=Wins_all_location_save, append=T, row.names=F, col.names=T,  sep=",")
}
############# SAVED ################


Wins_all_location_save <- "Wap_AllOcc_PlasmaEBIT_NoMixingCalc_CorrectedSubmission2022.csv"

#### Test
 test <- read_csv(Wins_all_location_save)
# test2 <- filter(test, Type_Good == FALSE)
#  
 
 
# # ########### PLOTTING STARTS HERE #########
# Wap0 <- read.csv(Wins_all_location_save)
# Wap <- Wap0
# #Wap <- filter(Wap0, B_WE != TRUE)
# 
# # GS
# Wap <- filter(Wap0, Type != "Large_L", !is.na(Rate_Pl_tot), !is.na(Rate_EBIT_tot), B_WE != TRUE)
# # Iso
# #Wap <- filter(Wap0, Type != "Large_L", !is.na(Rate_Pl_tot), !is.na(Rate_EBIT_tot), B_WE != TRUE,  Ei > 0)
# 
# Wap$Log_Rate_Pl_tot <- log10(Wap$Rate_Pl_tot)
# Wap$Log_Rate_EBIT_tot <- log10(Wap$Rate_EBIT_tot)

# #### HISTOGRAM
# ### try histogramming on the y axis
# p1 <- plot_ly(Wap, y=~Log_Rate_Pl_tot, type="histogram", name=~Type, color=~Type, colors="Set1") %>%
#    layout(
#      xaxis = list(title="Count"),
#      yaxis = list(title=FALSE, tick0=-10, dtick=2, showgrid=T, range=c(-18,10), showticklabels = FALSE)
#    )
# 
# ### PLASMA PLOT ### make horizontal with orientation='h'
# p2 <- plot_ly(Wap, x = ~M, y = ~Rate_Pl_tot, type = 'bar', name = ~Type, color=~Type, colors="Set1", showlegend=F) %>%
#   layout(yaxis = list(type="log", title = TeX("R_{\\alpha_{tot}} ... (Ion^{-1}s^{-1}))"), exponentformat='E',
#                        dtick=log10(1E2)),
#          xaxis = list(dtick =10, tickangle = -45, range=c(0, max(Wap$M)), title="Mass Number (A)")) %>%
#   config(mathjax = "cdn")
# 
# 
# ### 
# 
# 
# 
# 
# 
# fig <- subplot(p2, p1, widths = c(0.7, 0.3), titleY=TRUE, titleX=TRUE) %>%
#   layout(plot_bgcolor='#e5ecf6') %>%
#   config(mathjax = "cdn")
# fig

##### ...NEED TO MAKE shareY=T but add forced labelling of 1E(log10) on the right


### Set colour pallette
#cpal <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(8))

## Get rid of large multipolarity transitions
#Wap <- filter(Wap0, Type != "Large_L", Ei != 0, !is.na(Rate_Pl_tot), !is.na(Rate_EBIT_tot))
#Wap <- Wap0

# ### PLASMA PLOT ###
# plot_ly(Wap, x = ~M, y = ~Rate_Pl_tot, type = 'bar', name = ~Type, color=~Type, colors="YlOrRd") %>%
#   layout(yaxis = list(type="log", title = TeX("R_{\\alpha_{tot}} ... (Ion^{-1}s^{-1}))"), exponentformat='E'),
#          xaxis = list(dtick =10, tickangle = -45, range=c(0, 210), title="Mass Number (A)")) %>%
#   config(mathjax = "cdn")







# ### PLASMA PLOT ###
# plot_ly(Wap, x = ~M, y = ~Rate_EBIT_tot, type = 'bar', name = ~Type, color=~Type, colors="YlOrRd") %>%
#   layout(yaxis = list(type = "log", title = TeX("R_{EBIT} ... (Ion^{-1}s^{-1}))"), exponentformat='E'),
#          xaxis = list(dtick =10, tickangle = -45, range=c(0, 210), title="Mass Number (A)")) %>%
#   config(mathjax = "cdn")

#Wap_Bknown <- filter(Wap, !B_WE)

### try histogramming on the y axis
#h1 <- plotly(Wap)



######## WRITE RANKED DATA READY FOR TABULATION ######
# file1 <- "ranked_iso_plasma.csv"
# Wap <- arrange(Wap, desc(Rate_Pl_tot))
# write.table(Wap, file=file1, append=T, row.names=F, col.names=T,  sep=",")
# 
# file2 <- "ranked_iso_EBIT.csv"
# Wap <- arrange(Wap, desc(Rate_EBIT_tot))
# write.table(Wap, file=file2, append=T, row.names=F, col.names=T,  sep=",")


# threshold <- 1
# line.fmt = list(dash="dot", width = 1.5, color=NULL)
# fig_Nneectot <- fig_Nneectot %>% add_lines(y = threshold, line=line.fmt, name = "Observation Threshold", color="green")
# fig_Nneectot
# 
# 
# 
# add_trace(x = 0, y = 2/(pi*G), mode='text', text = TeX("\\frac{2}{\\pi \\Gamma_{neec}}"), textposition = 'left', textfont = list(size = 16)) %>%
#   add_trace(x = 0, y = 1/(pi*G), mode='text', text = TeX("\\frac{1}{\\pi \\Gamma_{neec}}"), textposition = 'left', textfont = list(size = 16)) %>%
#   add_trace(x = Eres, y = -0.04*1/(pi*G), mode='text', text = TeX("E_{res}"), textposition = 'below') %>%
#   add_segments(x = Eres-0.5*G, xend = Eres+0.5*G, y = 1/(pi*G), yend = 1/(pi*G), line=list(color = 'Black')) %>%
#   add_segments(x = Eres, xend = Eres, y = -0.02*1/(pi*G), yend = 2/(pi*G), line=list(color = 'Black', dash="dash")) %>%
#   add_segments(x = Eres, xend = 0, y = 2/(pi*G), yend = 2/(pi*G), line=list(color = 'Black', dash="dash")) %>%
#   add_segments(x = Eres-(0.5*G), xend = 0, y = 1/(pi*G), yend = 1/(pi*G), line=list(color = 'Black', dash="dash")) %>%
#   layout(yaxis = list(tickfont = list(color = "Blue"), title=TeX("Lor(E, \\Gamma_{neec})"), color = "Black", showticklabels=FALSE),
#          xaxis = list(tickfont = list(color = "Blue"), title = TeX("E (eV)"), showticklabels=FALSE), #categoryorder = "array", categoryarray = ~Irradiance_Wcm2um2),
#          showlegend = FALSE
#          # margin = list(
#          #   l = 50,
#          #   r = 50,
#          #   b = 100,
#          #   t = 100,
#          #   pad = 4
#          # )
#   ) %>%
#   config(mathjax = "cdn")









# fig_Rpl <- plot_ly(Wap, x = ~M, y = ~Rate_Pl_tot, type = 'bar', name = ~Type)
# fig_Rpl <- fig_Rpl %>% layout(yaxis = list(type = "log", title = "Plasma NEEC Rate s^-1"), xaxis = list(tickangle = -45))
# #fig_Nneectot
# # Add dotted threshold line
# #threshold <- 1
# #line.fmt = list(dash="dot", width = 1.5, color=NULL)
# #fig_Rpl <- fig_Rpl %>% add_lines(y = threshold, line=line.fmt, name = "Observation Threshold")
# fig_Rpl






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
  
