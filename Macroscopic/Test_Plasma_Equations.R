library(ggplot2)
source("Wins.R")
Mo <- read_csv("Mo_Additions_Rates2.csv")

#### Choose which facility
f <- 18
Te_eV <- 1.5E3
ne <- 1E24 ### e cm^-3
#Te <- 1.441
#<-22

### To get the average charge state
#Te_eV <- Te * 1E3
Z_name <- sprintf("Dependencies/CSD/CSD_Ave_Data/Z%003.f_AVE_CS.csv", 42)
Z_dat <- read_csv(Z_name)
Te_ind <- which.min(abs(FLYCHK_Te - Te_eV))[1]
Te_choice <- FLYCHK_Te[Te_ind]
ne_val <- ne
ne_ind <- which.min(abs(FLYCHK_ne - ne_val))[1]
ne_choice <- FLYCHK_ne[ne_ind]
CSD_File <- sprintf("Dependencies/CSD/CSD_fq_Data/CSD_%d.csv", 42)
CSD <- read_csv(CSD_File)
CSD_neTe <- filter(CSD, ne_cm3 == ne_choice, Te_eV == Te_choice)
Z_ave <- filter(Z_dat, Te == Te_choice, ne == ne_choice)$CS

#mass of ion # ... eV
mi_eV_c2 <- (Mo$Z[1]*mp + ( Mo$N[1] *mn ) + (Z_ave *me))

R_foc_m <- (facility$Spot_Diameter_um[f]*1E-6)/2

#depth of the plasma
dp_m <- c * facility$Pulse_Duration_fs[f] * 1E-15

tp_1 <- R_foc_m * sqrt(mi_eV_c2  / (Te_eV * Z_ave) ) * 1/c
tp_2 <- dp_m * sqrt(mi_eV_c2  / (Te_eV * Z_ave) ) * 1/c

facility <- mutate(facility, tpmax_1=double(1), tpmax_2=double(1) )



Z_name <- sprintf("Dependencies/CSD/CSD_Ave_Data/Z%003.f_AVE_CS.csv", 42)
Z_dat <- read_csv(Z_name)
CSD_File <- sprintf("Dependencies/CSD/CSD_fq_Data/CSD_%d.csv", 42)
CSD <- read_csv(CSD_File)

for(i in 1:length(facility$Facility)){
  ## Average charge state
  Te_eV <- facility$Te_chosen_keV[i] * 1E3
  Te_ind <- which.min(abs(FLYCHK_Te - Te_eV))[1]
  Te_choice <- FLYCHK_Te[Te_ind]
  ne_val <- facility$ne[i]
  ne_ind <- which.min(abs(FLYCHK_ne - ne_val))[1]
  ne_choice <- FLYCHK_ne[ne_ind]
  CSD_neTe <- filter(CSD, ne_cm3 == ne_choice, Te_eV == Te_choice)
  Z_ave <- filter(Z_dat, Te == Te_choice, ne == ne_choice)$CS
  
  #mass of ion # ... eV
  mi_eV_c2 <- (Mo$Z[1]*mp + ( Mo$N[1] *mn ) + (Z_ave *me))
  
  #focal length for tp1
  R_foc_m <- (facility$Spot_Diameter_um[f]*1E-6)/2
  #depth of the plasma for tp2
  dp_m <- c * facility$Pulse_Duration_fs[f] * 1E-15

  ## Assign data
  facility$tpmax_1[i] <-  R_foc_m * sqrt(mi_eV_c2  / (Te_eV * Z_ave) ) *1/c
  facility$tpmax_2[i] <- dp_m * sqrt(mi_eV_c2  / (Te_eV * Z_ave) ) * 1/c
    
}


##### Plot the plasma lifetime against irradiance
facility2 <- filter(facility, Irradiance_Wcm2um2 > 1E18)

ggplot(facility2, aes(x=Irradiance_Wcm2um2, y=tpmax_1, color=Facility)) +
  #geom_point() +
  geom_line()+
  #scale_x_log10() +
  #scale_y_log10()
  #geom_point(facility2, mapping=aes(x=Irradiance_Wcm2um2, y=tpmax_2))+
  geom_line()



