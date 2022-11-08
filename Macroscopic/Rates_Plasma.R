##########################
FLYCHK_ne <- c(1E12, 1E13, 1E14, 1E15, 1E16, 1E17, 1E18, 1E19, 1E20, 1E21, 1E22, 1E23, 1E24)
FLYCHK_Te <- c(0.5, 1, 1.5, 2, 5, 7, 10, 15, 23, 32, 52, 74, 100, 165, 235, 310, 390, 475, 655, 845, 1000, 1441, 1925, 2454, 3030, 3655, 4331, 5060, 5844, 6685, 7585, 8546, 10000, 20000, 50000, 100000)
#source("Macroscopic/Rates_Plasma.R")
#facility <- arrange(facility, Te_chosen_keV)

#source("Wins.R")
facility <- read_csv("Dependencies/Laser_Facilities/Big_Facilities_Irradiances_TeCalcd.csv")

MB <- function(E, Te){
  F_E <- 2 * ((E/(pi))^(1/2)) * ((1/(Te))^(3/2)) * exp(-E / Te)
  return(F_E)
}

#### All CSD data
fq_all <- list.files(path = "Dependencies/CSD/CSD_fq_Data",     # Identify all csv files in folder
                     pattern = "(CSD_\\d+.csv)", full.names = TRUE) %>%
  lapply(read_csv) %>%                                            # Store all files in list
  bind_rows             # Combine data sets into one data set



#####################################################################################
#####################################################################################
#################### PLASMA RATE   ##################################################
#####################################################################################
##### JUST SETTING Irradiance externally now via facility 
#beV, keV, eV
compute_Plasma_Rate <- function(Wins_all, FACIL_bool, facility_in, tau, INPUT_Te_Type, Te_input, ne_input, Ep_J, Rfoc_um, tlas_fs, Wavelen_nm, OUTPUT, OUTPUT2) {
  ### Microscopic parameters:
  #  i<-1000
  # Wins_all<-Wap[i,]; facility_in<- facility[j,]; INPUT_Te<-"Te"; Te_input<-1000*facility$Te_Beg_keV[j]; OUTPUT<-"N_neec_pw"; OUTPUT2<-"tot"
  ###### LASER PARAMETERS ######
  # first define the parameters we need for N_exc per pulse
  
  if(FACIL_bool == TRUE){
    #MACROSCOPIC PARAMETERS
    Epulse_J <- facility_in$Pulse_Energy_J
    Rfocal_m <- facility_in$Spot_Diameter_um*0.5E-6
    tlaser_s <- facility_in$Pulse_Duration_fs*1E-15
  } else if(FACIL_bool == FALSE) {
    Epulse_J <- Ep_J
    Rfocal_m <- Rfoc_um*1E-6
    tlaser_s <- tlas_fs*1E-15
  }
 
  Epulse_eV <- Epulse_J / e
  #median charge state .... need to use saha equation
  
  ###### NEEC WIDTH #####
  # Ar_1 <- compute_Ar("M1", A, Etrans_keV, 0.0078)
  # Ar_2 <- 0   ## Atomic rate, should be zero for this GS-GSA
  # Ar <- Ar_1 + Ar_2
  #Gamma_neec_J <- log(2) * (hbar / Wins_all$Thf_numeric)
  #Gamma_neec_J <- hbar * (Wins_all$Ar*(1+Wins_all$IC_tot))
  Gamma_neec_J <- hbar * (Wins_all$Ar*(1+Wins_all$ICC_sum))
  Gamma_d <- Gamma_neec_J / e
  Gamma_at <- 0
  Gamma_neec_eV <- Gamma_d + Gamma_at
  Lor_max <- 2/(pi*Gamma_neec_eV)
  
  
  #plasma lifetime
  #tp <- Rfocal_m * sqrt(mi_eV_c2 / (Te_eV * Z_ave))
  #plasma lifetime (mustnt forget to get rid of the c)
  #tp <- Rfocal_m *(1/c)*sqrt(mi_eV_c2 / (Te_eV * Z_ave)) #s
  tp <- tau # s (fix at 100ps)
  
  # #Electron temp from scaling law
  # if( INPUT_Te == "Scaling" ){
  #   Te_eV <-  facility$Te_Pond_keV*1000
  #   DELTA <- facility$Te_Pond_keV * 1000 * 4
  #   E_frac <- Wins_all$Gamma_neec_eV / DELTA
  # }
  ## If INPUT_Te_Type == "Te" this is ASTROPHYSICAL INPUT
  if( INPUT_Te_Type == "Te" ){
    #### OK WHATS THE PROPER FUCKING INPUT RANGE
    ### Ee and Te are always initially in keV
    Te_eV <- Te_input
    ### Just assume there are Z possible resonances for now
    E_frac  <-  MB((Wins_all$Ee*1000), Te_eV) * (Gamma_neec_eV)
    MB_val <- MB((Wins_all$Ee*1000), Te_eV)
    ne <- ne_input
    ne_actual <- ne
    Vp <- 1
    
    
  } else if( INPUT_Te_Type == "Params" ) {
    # Epulse_J <- Ep_J
    # Rfocal_m <- Rfoc_um*1E-6
    # tlaser_s <- tlas_fs*1E-15
    P_W <-  Epulse_J/(tlaser_s)
    I_Wcm2 <- P_W / (pi * (1E2*Rfocal_m)^2)
    Wavelen_um <- Wavelen_nm*1E-3
    #wavelength_um <- facility$Wavelength_nm[i]*1E-3
    Irr_Wcm2um2 <- ((Wavelen_nm*1E-3)^2) * I_Wcm2
    ### Laser intensities for scaling laws
    I16 <- I_Wcm2/1E16
    I18 <- I_Wcm2/1E18
    
    #Short scale length profiles 
    SSLP_keV <- 8*(I16 * Wavelen_um^2)^(1/3)
    #beg scaling
    Beg_keV <-  230*(I18 * Wavelen_um^2)^(1/3)
    #ponderomotive scaling
    Pond_keV <- 3.6* I16 * Wavelen_um^2
    
    if(Irr_Wcm2um2 < 3.31E16){
      Te_chosen_keV <- Pond_keV
    } else if(Irr_Wcm2um2 < 5.1E17){
      Te_chosen_keV <- Pond_keV
    } else {
      Te_chosen_keV <- Beg_keV
    }
    
    Te_eV <- Te_chosen_keV*1000
    E_frac  <-  MB((Wins_all$Ee*1000), Te_eV) * (Gamma_neec_eV)
    MB_val <- MB((Wins_all$Ee*1000), Te_eV)
    f <- 0.5
    
    ##to get elec number density
    Ne <- f*Epulse_eV / Te_eV 
    dp <- c*100*tlaser_s #cm
    Vp <- pi*((1E2*Rfocal_m)^2)*dp #cm^3
    ne <- Ne/Vp ## e's cm^-3
    ne_actual <- ne
  
  }
  # } else if( INPUT_Te == "Astro" ){
  #    Te_chosen_keV <- 
  # }
  # 
  
  
  ##### ABSORBTION FRACTION #####
  # for(i in 1:length(facility$tp)){
  # facility$f[i] =(1.2E-15) * facility$Intensity_Wcm2[i] * 10^(0.74)
  # }
  # if(f>0.5){
  #   f <- 0.5
  #}
#### MAY NEED UNCOMMENTING ####  
  # f <- 0.5
  # 
  # ##to get elec number density
  # Ne <- f*Epulse_eV / Te_eV
  # dp <- c*100*tlaser_s #cm
  # Vp <- pi*((1E2*Rfocal_m)^2)*dp #cm^3
  # ne <- Ne/Vp ## e's cm^-3
  # ne_actual <- ne
  
  
  #MUST REPLACE WITH Wins_all$Z
  #Wins_allZ <- 40
  Te_eV <- Te_eV
  
  
  ##### Read the correct CSD fq file ####
  
  
  if((Wins_all$Z < 79) && ((1E-3*Te_input) <= 100)){
    #### READ THE CORRECT Z_ave ####
    Z_name <- sprintf("Dependencies/CSD/CSD_Ave_Data/Z%003.f_AVE_CS.csv", Wins_all$Z)
    
    Z_dat <- read_csv(Z_name)
    
    ##### Read the correct CSD fq file ####
    
    
    ## Round to nearest ne and Te
    #### test choose nearest val from Te and ne veactors
    # FLYCHK_ne <- c(1E12, 1E13, 1E14, 1E15, 1E16, 1E17, 1E18, 1E19, 1E20, 1E21, 1E22, 1E23, 1E24)
    # FLYCHK_Te <- c(0.5, 1, 1.5, 2, 5, 7, 10, 15, 23, 32, 52, 74, 100, 165, 235, 310, 390, 475, 655, 845, 1000, 1441, 1925, 2454, 3030, 3655, 4331, 5060, 5844, 6685, 7585, 8546, 10000, 20000, 50000, 100000)
    # 
    #Te_eV <- 880
    
    #Te
    Te_ind <- which.min(abs(FLYCHK_Te - Te_eV))[1]
    Te_choice <- FLYCHK_Te[Te_ind]
    #Te_eV <- Te_choice
    
    
    ne_val <- ne
    ne_ind <- which.min(abs(FLYCHK_ne - ne_val))[1]
    ne_choice <- FLYCHK_ne[ne_ind]
  
    #CSD_File <- sprintf("Dependencies/CSD/CSD_fq_Data/CSD_%d.csv", Wins_all$Z)
    
    #CSD <- read_csv(CSD_File)
    CSD <- filter(fq_all, Z==Wins_all$Z)
    CSD_neTe <- filter(CSD, ne_cm3 == ne_choice, Te_eV == Te_choice)
  
    #### We read the charge state of the CS chosen by the closest Te and ne selected by the laser params
    Z_ave <- filter(Z_dat, Te == Te_choice, ne == ne_choice)$CS
    fq <- CSD_neTe$frac[which(CSD_neTe$q==(Wins_all$CS+1))]
    #if(Wins_all$Z >)
    
    
  } else {
    Z_ave <- Wins_all$Z
    fq <- 1
  } 
  
  #mass of ion # ... eV
  # mi_eV_c2 <- (Mo$Z*mp + ( Mo$N *mn ) + (Z_ave *me))
  
  ### SPECIES 
  #Ee <- 1000 #keV
  #vres <- sqrt(2*(Ee*1E3)/(me))  #c
  #vres <- 100*c*vres   #cm s^-1
  
  vres_rel <- sqrt(1 - (1/((1+(Wins_all$Ee*1E3) / me))^2))  #c
  vres <- 100*c*vres_rel #cm s^-1
  #vres <- sqrt(2*(Wins_all_plasma$Ee[1]*1E3)/(me)) #c
  
  ## resonant flux
  phi_res <- ne_actual * vres #cm^-2 s^-1
  
  #Total ICC  
  if(OUTPUT2 == "tot"){
    Sexp <- Wins_all$Stot * 1E-24  #cm^2 eV
  }
  #Shell ICC
  if(OUTPUT2 == "p"){
    #Sexp <- Wins_all$S * Wins_all$frac_CS * 1E-24  #cm^2 eV
    Sexp <-  Wins_all$S * 1E-24 * fq   #cm^2 eV  #* Wins_all$Gamma_Scale_Factor
  } 
  

  ## resonant cross section
  sigma_res <- Sexp*Lor_max 
  
  ## Macroscopic NEEC RATE [s^-1] based on equiprobable electron number density
  ###### NEEDED MULTIPLYING BY THE RESONANT FRACTION
  ni <- (ne_actual / Z_ave) 
  Ni <- ni * Vp
  #R_neec <- phi_res * sigma_res  *  Ni * E_frac
  R_neec <- phi_res * MB_val * Sexp * (2/pi) #* Ni
  
  
  ### Lifetime of the plasma and CSD in the plasma
  
  
  ############ CHECKING AGAINST GUNST  ###################
  # Te_eV <-  facility$Te_Beg_keV*1000
  # Z_ave
  # mi_eV_c2 <- (Z*mp + ( (M - Z) *mn ) + (Z_ave *me))
  # #mi_Kg <- mi_eV_c2 * e / (c^2)
  # #plasma lifetime
  # tp <- Rfocal_m * sqrt(mi_eV_c2 / (Te_eV * Z_ave))  
  
  ############ CHECKING AGAINST GUNST  ###################
  
  
  ## NEEC's per pulse
  N_neec_pp <- R_neec * tp * Ni
  
  ## NEEC's per week
  N_neec_pw <- N_neec_pp * (facility_in$Max_Repetition_Rate_Hz * 60*60*24*7)
  
  if(OUTPUT == "fq"){
    return(fq)
  }
  if(OUTPUT == "N_neec_pw"){
    return(N_neec_pw)
  }
  if(OUTPUT == "Rate"){
    return(R_neec)
  }
  if(OUTPUT == "Ne"){
    return(Ne)
  }
  if(OUTPUT == "Te"){
    return(Te)
  }
  if(OUTPUT == "Nion"){
    return(Ni)
  }
  
}

#compute_Plasma_Rate(Wins_all_test[139,], facility[9,])

##########################      saha <- function()
##### Look at the degeneracy of the levels for the sequential value based on a spectrum form a charge state


