library(tidyverse)
library(stringr)

source("Wins.R")

Sys.setenv("plotly_username" = "bw979")
Sys.setenv("plotly_api_key" = "v7KDKKj5jCkkpKsid0kS")

# #### Load the 84Rb resonance spectrum from the master NEEC data
# n_index <- read_csv("App/n_to_ICshell_Conversion.csv")
# atomic_colours <- read_csv("App/Atomic_Colours.csv")
# wap0 <- read_csv("App/All_RatesCalcd_Viva.csv")
# Mo <- read_csv("App/Mo_Additions_Rates2.csv")
# wap <- bind_rows(wap0, Mo)
# wap <- select(wap, -Rneec_Pl, -Rneec_SEBIT) %>%
#   filter(Z>9)
# wap_N0 <- filter(wap, Ei==0)
# wap_Ni <- filter(wap, Ei!=0)

### 84Rb isomer transition
#Cand <- read_csv("Rb_Addition.csv")
### 93Mo case
Cand <- read_csv("Mo_Additions_Rates2.csv")
#Cand$S2 <- Cand$S/Cand$Gamma_Scale_Factor
##Use unscaled ICC's
##### NO NEED FOR THIS AS RECALCULATING ICC AND S LATER ANYWAY
# for(i in 1:length(Cand$S)){
#   #Type, ss,  Zin, Q, n, OUTPUT, MIX, MR, MIX_type)
#   ### Recalculate unscaled S (in wap it has been scaled by the gamma detection factor)
#   Cand$ICC2[i] <- Compute_ICC_singles(Cand$Type[i], Cand$subshell[i], Cand$Z[i], Cand$Q[i], Cand$Occ[i], "p", FALSE, NA, NA)
#   Cand$S2[i] <- compute_S(8.5, Cand$Jat[i], 10.5, Cand$ICC2[i], Cand$Ar[10], Cand$Q[i], Cand$Ebind[i])
# }


Cand <- select(Cand, -Rate_Pl, -Rate_SEBIT, -Npw_Pl_S, - Npw_SEBIT)

### Read Energy loss array ... ENERGIES ARE ALWAYS IN MeV and MeV/um
#dedx <- read_csv("Dependencies/Beam_Data/E_dEdx.csv")
dedx <- read_csv("Dependencies/Beam_Data/Mo_C_dEdx.csv")

# Cand$Ji <- "6-"
# Cand$Ji_double <- 6
# Cand$Jf <- "5-"
# Cand$Jf_double <- 5
# Cand$Ei <- 463.59
# Cand$Ef <- 466.64
# Cand$Btype <- "BM1W"
# Cand$B <- 0.08
# Cand$Type <- "M1"
# Cand$B2type <- "BE2W"
# Cand$B2 <-
# Cand$Thi <- "20.26 min"
# Cand$Thi_numeric <- 1215.6
# Cand$Thf <- "9ns"
# Cand$Thf_numeric <- 9E-9
# ## Set transition Energy
# Etrans <- 4
# Cand$Q <- Etrans
Cand$Ee <- Cand$Q - Cand$Ebind
Cand$Ebeam <- mp.me * Cand$Ee * 1E-3 * 84   # MeV

#Effective resonance energy
Ebeam_eff <- Cand$Q[1] - mean(Cand$Ebind)

## Set target thickness
Thickness <- 100 #um
stepsize <- 0.001
nsteps <- (Thickness/stepsize) + 1
## Set Entrance beam energy
#Ebeam <- 7.35 #MeV/u 
#Ebeam_MeV <- 84*Ebeam

## Entrance Beam Energy
Ebeam_MeV <- max(Cand$Ebeam) #MeV

### bulk values
Ib <- 1E5       #Beam intensity per second
#f_ch <- 0.001     #channeled fraction

Z <- 42
Mc <- 12.011 # g/mol
Na <- 6.02E23 #Avagadros number atoms/mol
rho_c <- 2.26 # gcm-3
Nc <- rho_c * (1/Mc) * Na
ne <- Z* Nc  #e- cm^-3

#n_e <- 1E23  #e- cm^-3
ne <- ne * 1E-12  #e- um^-3
n_e<-ne


### Interpolate dEdx array
Cand <- mutate(Cand, dedx_res = double(1), P_neec = double(1), fq = 1, R_neec = double(1), x = double(1))

### Resonance Strength Calculation
for(i in 1:length(Cand$AX)){
  #Wap <- Cand
  ##### CALCULATE Ar, ICC and Stot
#  if(  !is.na(Wap$B[i]) && !is.na(Wap$Ji_double[i]) && !is.na(Wap$Jf_double[i]) && !is.na(Wap$Type[i]) && Wap$Type[i] != "E0" && Wap$Type[i] != "Large_L") {
    # ######Compute_Ar... reuturns in SI [s^-1] ... also trying ICC_tot
    # #if its mixed
    # if(  (Wap$MR[i]!=0) && (!is.na(Wap$MR[i])) ){
    #   Cand$Ar[i] <- compute_Ar(Wap$Type[i], Wap$M[i], Wap$Q[i], Wap$B[i], TRUE, Wap$MR[i] )
    #   ## If its "M1 + E2" mixed
    #   if( !is.na(Wap$MIX_type[i]) && Wap$MIX_type[i] == "M1+E2"  ){
    #     try(Cand$ICC_tot[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 1, "tot", TRUE,  Wap$MR[i], Wap$MIX_type[i])  )
    #     try(Cand$ICC[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "p", TRUE, Wap$MR[i], Wap$MIX_type[i])  )
    #     try(Cand$ICC_sum[i] <- Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], Wap$Occ[i], "SUMp", TRUE, Wap$MR[i], Wap$MIX_type[i])  )
    #   }
    # } else {
      Cand$Ar[i] <- compute_Ar(Cand$Type[i], Cand$M[i], Cand$Q[i], Cand$B[i], FALSE, Cand$MR[i] )
      #Compute_ICC_singles(Cand$Type[i], Cand$subshell[i], Cand$Z[i], Cand$Q[i], Cand$Occ[i], "p", FALSE, NA, NA)
      try(Cand$ICC_tot[i] <- Compute_ICC_singles(Cand$Type[i], Cand$subshell[i], Cand$Z[i], Cand$Q[i], Cand$Occ[i], "tot", FALSE, NA, NA) )
      try(Cand$ICC[i] <- Compute_ICC_singles(Cand$Type[i], Cand$subshell[i], Cand$Z[i], Cand$Q[i], Cand$Occ[i], "p", FALSE, NA, NA) )
      #try(Cand$ICC_sum[i] <- Compute_ICC_singles(Cand$Type[i], Cand$Z[i], Cand$Q[i], Cand$Occ[i], "SUMp", FALSE, 0, "") )
#    }
#  }
  
  
  #   # #Compute_S Jd, Jat,Ji,ICC, Ar, Etrans_keV, Ebind_keV
  #   # Wap$S[i]  <-  compute_S(Wap$Jf_double[i], Wap$Jat[i], Wap$Ji_double[i], Wap$ICC[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
  #   # Wap$S_tot[i] <- compute_S(Wap$Jf_double[i] , Wap$Jat[i], Wap$Ji_double[i], Wap$ICC_tot[i], Wap$Ar[i], Wap$Q[i], Wap$Ebind[i])
  # 
  # 
  # ### Total depletion level width
  Cand$Gamma_Thalf_eV[i] <- (hbar/e) * (log(2) / Cand$Thf_numeric[i])
  Cand$Gamma_eV[i] <-  (hbar/e) * (Cand$Ar[i]*(1+Cand$ICC_sum[i]))
  Cand$Gamma_Scale_Factor[i] <- 1/(1+Cand$ICC_sum[i])

  ### COMPUTE THE RESONANCE STRENGTH
  print(i)
  Cand$Stot[i]  <-  compute_S(Cand$Jf_double[i], Cand$Jat[i], Cand$Ji_double[i], Cand$ICC_tot[i], Cand$Ar[i], Cand$Q[i], mean(Cand$Ebind))#*Cand$Z[i] 
  try(Cand$S[i]  <-  compute_S(Cand$Jf_double[i], Cand$Jat[i], Cand$Ji_double[i], Cand$ICC[i], Cand$Ar[i], Cand$Q[i], Cand$Ebind[i]) ) #*Cand$Gamma_Scale_Factor[i] )
  #try(Cand$S_sum[i]  <-  compute_S(Cand$Jf_double[i], Cand$Jat[i], Cand$Ji_double[i], Cand$ICC_sum[i], Cand$Ar[i], Cand$Q[i], Cand$Ebind[i]) )
  
}



#### PROBABILITY AND NEEC RATE
for(i in 1:length(Cand$AX)){
  print(i)
  E_val <- Cand$Ebeam[i] # MeV
  E_ind <- which.min(abs(dedx$E_Ion - E_val))[1]
  dedx_choice <- dedx$dEdXTot[E_ind]*1E-3  ## MeV/um
  Cand$dedx_res[i] <- dedx_choice
  Cand$P_neec[i] <- ne * Cand$fq[i] * Cand$S[i] *1E-24 * 1E8 * 1E-6 *Cand$M[i]* mp.me * 1/(Cand$dedx_res[i])
  Cand$R_neec[i] <- Cand$P_neec[i] * Ib
}

P_NEEC <- sum(filter(Cand, !is.na(P_neec))$P_neec)

#### DOING THE alphaTot single calculation
## Finds out nearest dedx known
E_val <- Ebeam_eff # MeV
E_ind <- which.min(abs(dedx$E_Ion - E_val))[1]
dedx_choice <- dedx$dEdXTot[E_ind]*1E-3  ## MeV/um
dedx_res_alphaTot <- dedx_choice

P_NEEC_alphaTot <- ne * compute_S(Cand$Jf_double[1], mean(Cand$Jat), Cand$Ji_double[1], Cand$ICC_tot[1], Cand$Ar[1], Cand$Q[1], mean(Cand$Ebind)) *1E-24 * 1E8 * 1E-6 *Cand$M[1]* mp.me * 1/(dedx_res_alphaTot)


### SUM Rates and PLOT
#Cand_2 <- filter(Cand, !is.na(R_neec), shell != "L")
#rneec <- signif(sum(Cand_2$R_neec), 5)

##### Array of energies explored by the beam due to this thickness
E_explore <- tibble(x = seq(0, Thickness, stepsize), Eb = double(nsteps), S=double(1), shell = double(1))
for(i in 1:length(E_explore$x)){
  #i<-2
  if(i==1){
    E_explore$Eb[i] <- Ebeam_MeV
  } else {
    E_val <- E_explore$Eb[i-1] # MeV
    E_ind <- which.min(abs(dedx$E_Ion - E_val))[1]
    dedx_choice <- dedx$dEdXTot[E_ind]*1E-3  ## MeV/um
    E_explore$Eb[i] <- E_explore$Eb[i-1] - (dedx_choice * stepsize) 
  }
  
 # E_explore$S <- Cand$S[which(Cand$x == E_explore$x[i])]
  
}    

### How far through the target are these resonances?
for(i in 1:length(Cand$AX)){
  E_val <- Cand$Ebeam[i] #MeV
  E_ind <- which.min(abs(E_explore$Eb - E_val))[1] #MeV
  x_choice <- E_explore$x[E_ind]  # um
  Cand$x[i] <- x_choice
}


#S_vals = c(rep(0, length(E_explore$x)))


for(i in 1:length(E_explore$x)){
  #i <-1
  if(E_explore$x[i] %in% Cand$x){
    E_explore$S[i] <- Cand$S[which(Cand$x == E_explore$x[i])[1]]
    E_explore$shell[i] <- Cand$shell[which(Cand$x == E_explore$x[i])[1]]
    #S_vals[i] <- Cand$S[which(Cand$x == E_explore$x[i])][1]
    #E_explore$E_fac <- E_explore$
  } else {
    E_explore$S[i] <- 0
    E_explore$shell <- NA
  }
}

#E_explore$E_fac <- E_explore$
  


# ### SUM Rates and PLOT
# Cand_2 <- filter(Cand, !is.na(R_neec)) #shell != "L")
# rneec <- signif(sum(Cand_2$R_neec), 5)

## Detection
f_det <- 0.1


# ### Resonance Spectrum
# p1 <- plot_ly(Cand, x = ~Ebeam, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% #, marker = list(color = label_colours)) %>%
#   layout(yaxis = list(type = "log", exponentformat = "E", title = "NEEC Resonance Strength / beV"), #range = c(1E-10, max(wap_filtered$S))),
#          xaxis = list(tickangle = -45, title="Resonance Beam Energy (MeV/u)", range = c(min(Cand_2$Ebeam), max(Cand_2$Ebeam))),
#          legend=list(name=~shell, showlegend=TRUE)
#          # yaxis2 = list(
#          #   tickfont = list(color = "red"),
#          #   overlaying = "y",
#          #   side = "right",
#          #   title = "MB" ,
#          #   showlegend = F)
#   )
# 
# p1
# 
# 
# #### NEW ARRAY BE_S
# # step <- 0.1
# # steps <- (Thickness/step) + 1
# # BE_S <- tibble(x = seq(0, Thickness, stepsize), Eb = double(steps))
# 
# ### Resonance depth
# plot_ly(Cand, x = ~x, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% #, marker = list(color = label_colours)) %>%
#   #add_trace(E_explore, x = ~x, y = ~, type = 'bar', type = 'scatter') %>%
#   layout(yaxis = list(type = "log", exponentformat = "E", title = "NEEC Resonance Strength / beV"), #range = c(1E-10, max(wap_filtered$S))),
#          xaxis = list(tickangle = -45, title="Target Depth (um)", range = c(min(E_explore$x), max(E_explore$x))),
#          legend=list(name=~shell, showlegend=TRUE) 
#          # yaxis2 = list(
#          #   tickfont = list(color = "red"),
#          #   overlaying = "y",
#          #   side = "right",
#          #   title = "MB" ,
#          #   showlegend = F)
#          
#   )#    %>%  

### Plot beam energy 
### Resonance depth
#E_explore <- filter(is.even(E_explore$x))


# # getting number of rows in R
# rows <- nrow(E_explore)
# # extracting odd rows 
# odd_rows <- seq_len(rows) %% 4
# # getting data from odd data frame
# data_mod <- E_explore[odd_rows == 1, ]

# p2 <- plot_ly(E_explore, x = ~x, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) #%>% #, marker = list(color = label_colours)) %>%
#   #add_trace(E_explore, x = ~x, y = ~Eb, type = 'scatter', mode='lines') %>%
#   layout(yaxis = list(type = "log", exponentformat = "E", title = "Beam Energy (MeV)"), #range = c(1E-10, max(wap_filtered$S))),
#          xaxis = list(tickangle = -45, title="Target Depth (um)", range = c(min(E_explore$x), max(E_explore$x))),
#          legend=list(name=~shell, showlegend=TRUE)
#         # legend=list(name=~shell, showlegend=TRUE) 
#          # yaxis2 = list(
#          #   tickfont = list(color = "red"),
#          #   overlaying = "y",
#          #   side = "right",
#          #   title = "MB" ,
#          #   showlegend = F)
#          
#   )    %>%
#   add_lines(x = ~x, y = ~Eb, type = 'scatter', mode = 'lines', color="#DE3163", showlegend=T, name="Beam Energy") 
# 
# p2

# #### DOUBLE AXES IRRADIANCE AND
# y2 <- list(
#   tickfont = list(color = "red"),
#   overlaying = "y",
#   side = "right",
#   title = "second y axis",
#   type = 'log',
#   color = 'black'
# )

#################################################
####### CURRENTLY THESIS STANDARD PLOT ##########
#################################################
##### ENERGY EXPLORED AND RESONANCE STRENGTH PLOT
plot_ly(Cand, x = ~x, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% #, marker = list(color = label_colours)) %>%
  layout(plot_bgcolor='#e5ecf6',
         yaxis = list(type = "log", exponentformat = "E", title = "NEEC Resonance Strength (beV)" , side="right"  ), #range = c(1E-10, max(wap_filtered$S))),
         yaxis2 = list(
           tickfont = list(color = "red"),
           overlaying = "y",
           side = "left",
           title = "Energy of Ion Inside Target (MeV)" ,
           showlegend = F,
           range = c(0, Ebeam_MeV)),
         xaxis = list(tickangle = -45, title="Target Depth (um)", range = c(0, 80), dtick = 10, showgrid=T),
         legend=list(name=~shell, showlegend=TRUE)
         ) %>%
  add_lines(x = ~E_explore$x, y = ~E_explore$Eb, type = 'scatter', mode = 'lines', color=list("black"), yaxis="y2", showlegend=T, name="Ion Energy ") %>%
  config(mathjax = "cdn") %>%
  layout(legend= list(x=0.8, y=0.9))









#plotly_IMAGE(p2, format = "png", out_file = "Beam_Energy")
#htmlwidgets::saveWidget(p2, "Beam_Energy.html")

#plot(E_explore$x, E_explore$Eb)

# ### E, depth
# 
# #### DOUBLE AXES Resonance Strength and Beam Energy
# y2 <- list(
#   tickfont = list(color = "red"),
#   overlaying = "y",
#   side = "right",
#   title = "second y axis",
#   type = 'log',
#   color = 'red'
# )
# 
# fig4a <- plot_ly(Cand_2, x = ~x, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% 
#   add_trace(y = ~ , yaxis = 'y2', type = 'scatter')
# fig4a <- fig4a %>% layout(yaxis = list(tickfont = list(color = "blue"), type = 'log', title="Laser Repetition Rate, Hz", color = 'blue' ),
#                           yaxis2 = list(tickfont = list(color = "orange"), overlaying = "y", side = "right", title = "Irradiance, Wcm^-2 um^2", type = 'log', color = 'orange'),
#                           xaxis = list(title = "Facility", categoryorder = "array", categoryarray = ~Irradiance_Wcm2um2),
#                           showlegend = FALSE, 
#                           pad = 100,
#                           margin = list(
#                             l = 50,
#                             r = 50,
#                             b = 100,
#                             t = 100,
#                             pad = 4
#                           )
# ) 
# 
# fig4a


##### Array of energies explored by the beam due to this thickness


# ### MAX Thickness required
# xmax <- max(Cand$x)
# 
# ##### Coulex Rate
# sig_coul <- 0.1E-28    #[m^2] ... = 100mb
# rho_nuc <- 5E28    #m^-3
# deltaZtot <- Thickness * 1E-6
# rcoul <- Ib*sig_coul*rho_nuc*deltaZtot
# 
# p <- sprintf("Rneec: %.3g       Rcoul: %.3g         Max Thickness: %.3g  um",rneec,rcoul,xmax)
# print(p)