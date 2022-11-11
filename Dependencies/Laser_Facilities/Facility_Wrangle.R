#source("Wins.R")
facility <- read_csv("Dependencies/Laser_Facilities/All_Facilities_Irradiances.csv")

MB <- function(E, Te){
  F_E <- 2 * ((E/(pi))^(1/2)) * ((1/(Te))^(3/2)) * exp(-E / Te)
  return(F_E)
}

for(i in 1:length(facility$Power_W)){
  #i<-9
  Epulse_J <- facility$Pulse_Energy_J[i]
  Rfocal_m <- facility$Spot_Diameter_um[i]*0.5E-6
  tlaser_s <- facility$Pulse_Duration_fs[i]*1E-15
  
  facility$Power_W[i] <-  Epulse_J/(tlaser_s)
  facility$Intensity_Wcm2[i] <- facility$Power_W[i] / (pi * (1E2*Rfocal_m)^2)
  
  #wavelength_um <- facility$Wavelength_nm[i]*1E-3
  facility$Irradiance_Wcm2um2[i] <- ((facility$Wavelength_nm[i]*1E-3)^2) * facility$Intensity_Wcm2[i]
  
  ### Laser intensities for scaling laws
  I16 <- facility$Intensity_Wcm2[i]/1E16
  I18 <- facility$Intensity_Wcm2[i]/1E18
  wavelength_um <- facility$Wavelength_nm[i]*1E-3
  #Short scale length profiles 
  facility$Te_SSLP_keV[i] <- 8*(I16 * wavelength_um^2)^(1/3)
  #beg scaling
  facility$Te_Beg_keV[i] <-  230*(I18 * wavelength_um^2)^(1/3)
  #ponderomotive scaling
  facility$Te_Pond_keV[i] <- 3.6* I16 * wavelength_um^2
  
  if(facility$Irradiance_Wcm2um2[i] < 3.31E16){
    facility$Te_chosen_keV[i] <- facility$Te_Pond_keV[i] 
  } else if(facility$Irradiance_Wcm2um2[i] < 5.1E17){
    facility$Te_chosen_keV[i] <- facility$Te_Pond_keV[i]
  } else {
    facility$Te_chosen_keV[i] <- facility$Te_Beg_keV[i]
  }
  
  
  # ############ CHECKING AGAINST GUNST  ###################
   Z <- 42
   M <- 93
   f<-0.5
   Te_eV <-  facility$Te_chosen_keV[i]*1000
   Epulse_eV <- Epulse_J / e
   Z_ave <- Z/2
  mi_eV_c2 <- (Z*mp + ( (M - Z) *mn ) + (Z_ave * me))
  # #mi_Kg <- mi_eV_c2 * e / (c^2)
  # #plasma lifetime
  # 
  # 
  # 
  ##to get elec number density
  Ne <- f*Epulse_eV / Te_eV
  dp <- c*100*tlaser_s #cm
  Vp <- pi*((1E2*Rfocal_m)^2)*dp #cm^3
  ne <- Ne/Vp ## e's cm^-3
  
  facility$Ne[i] <- Ne
  ### THIS IS IN mc^-1 you idiot need to convert to s by diving by c
  tp <- Rfocal_m *(1/c)*sqrt(mi_eV_c2 / (Te_eV * Z_ave)) #s

  facility$tp[i] <- tp
  facility$ne[i] <- Ne/Vp
  ############ CHECKING AGAINST GUNST  ###################
  facility <- facility[order(facility$Te_chosen_keV),]
}
if(!file.exists("Dependencies/Laser_Facilities/Big_Facilities_Irradiances_TeCalcd.csv")){
  write.table(facility, file="Dependencies/Laser_Facilities/Big_Facilities_Irradiances_TeCalcd.csv", append=T, row.names=F, col.names=T,  sep=",")
}


