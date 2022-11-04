library(tidyverse)
library(ggplot2)
library(stats)
library(grid)
library(plotly)

source("Wins.R")
## No Mo isomer B values in this one:
#wap <- read_csv("Outputs/MASTER/All/All_RatesCalcd_Viva.csv")

### Read in data
#Wap_N0 <- read_csv("Outputs/MASTER/All/AnN0_all.csv")
#Wap_Ni <- read_csv("Outputs/MASTER/All/AnNi_all.csv")
#wap <- bind_rows(Wap_N0, Wap_Ni)
additions <- read_csv("Mo_Additions_Rates2.csv")
wap <- read_csv("Outputs/MASTER/All/All_RatesCalcd_Viva.csv")
wap <- bind_rows(wap, additions)

######## For 57Fe  ######################
wap_f <- filter(wap, AX == "57FE", Ef == 14.4129)

######## For 93Mo ########################
wap_f2 <- filter(wap, AX == "93MO", Q < 4.86)

######## For 93Mo ########################
wap_f3 <- filter(wap, AX == "73GE", Q==13.2845)

#wap_f2$B <- 3.5
## Read the Abinit data 
Abinit <- read.csv("Microscopic/Abinit_Apr2022.csv")
Abinit_GS <- filter(Abinit, GS==TRUE)


MB <- function(E, Te){
  F_E <- 2 * ((E/(pi))^(1/2)) * ((1/(Te))^(3/2)) * exp(-E / Te)
  return(F_E)
}

############# SET ENERGIES AND TEMPS FOR PLOTTING  ########
##Electron Temperature
Te <- 12000  #eV
##Number Density
ne <- 1E29  #cm^-3

## Impact energy range
#E <- seq(0.1,max(Abinit$Eres_eV*1E-3),0.1) #keV
E <- seq(0.01,max(wap_f3$Ee),0.01) #keV
v_rel <- sqrt(1 - (1/((1+(E*1E3) / me))^2)) #c
v <- 100*c*v_rel #cms^-1

### NEED TO EXPRESS h IN THE CORRECT UNITS ###
## Effective DOS
g <-  2*(2*pi*me_cms*Te/(h_eV^2))^(3/2)
g_adj <- g*(2/sqrt(pi))

## Reduced chemical potential
mu_Te <- ((3*sqrt(pi) * ne)/(4*g_adj))^(2/3)

mu <- mu_Te*Te

## Actual DOS
#ge <-((2*me_cms)^(3/2))*(1/(2*(pi^2)*(hbar_eV^3)))*sqrt(E*1E3)
ge <- ((4)/(sqrt(pi))) * (((2*pi*me_cms)/(h_eV^2))^(3/2))*((E*1E3)^0.5)

## Electron distribution functions
# MB
F_MB <- MB(E*1E3, Te)
# FD
F_FD <- 1/(1+((exp(E*1E3/Te)*exp(-mu_Te))))


## Maxwell Boltzmann flux
phi_MB <- ne * F_MB * v

## Fermi dirac Flux
phi_FD <- F_FD * ge * v




# ############ LOOK AT 57Fe #################
# ## Select transition
#wap_f <- filter(wap, AX == "57FE", Ef == 14.4129)
#
# ## Graph
# plot_ly(wap_f, x = ~Ee, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% #, marker = list(color = label_colours)) %>%
#   layout(yaxis = list(type = "log", exponentformat = "E", title = "57Fe NEEC Resonance Strength / beV"), #range = c(1E-10, max(wap_filtered$S))),
#          xaxis = list(tickangle = -45, title="Electron Impact Energy / keV", range = c(min(E), max(E))),
#          legend=list(name=~shell, showlegend=TRUE)
#          # yaxis2 = list(
#          #   tickfont = list(color = "red"),
#          #   overlaying = "y",
#          #   side = "right",
#          #   title = "MB" ,
#          #   showlegend = F)
# 
#   )  %>%
#   add_lines(x = ~E, y = ~phi_MB, type = 'scatter', mode = 'lines', showlegend=T, name="Phi_MB") %>%
#   add_lines(x = ~E, y = ~phi_FD, type = 'scatter', mode = 'lines', showlegend=T, name="Phi_FD")
# 

############# 93Mo ##################

Abinit_GS$Ee <- Abinit_GS$Eres_eV*1E-3

MBFD <- tibble(E=E, phi_FD=phi_FD, phi_MB=phi_MB)

## Add in a NEEC cross section for comparison with Palffy
#Abinit_GS$CS <-Abinit_GS$Spdb*

# ## Graph
# plot_ly(data=MBFD, x = ~E, y = ~phi_FD, type = 'scatter', mode="lines", showlegend=TRUE, name=TeX("\\phi_{FD}")) %>% #, marker = list(color = label_colours)) %>%
#   layout(yaxis = list(type='log', exponentformat = "E", title = TeX("\\phi(E) cm^{-2}s^{-1}eV^{-1}")), #range = c(1E-10, max(wap_filtered$S))),
#          xaxis = list(tickangle = -45, title="Electron Impact Energy  (keV)", range = c(min(E), max(E))),
#          legend=list(showlegend=TRUE),
#          yaxis2 = list(
#                       tickfont = list(color = "red"),
#                       overlaying = "y",
#                       side = "right",
#                       title = "Energy of Ion Inside Target (MeV)" ,
#                       showlegend = F)#,
#                       #range = c(0, Ebeam_MeV))
#                       
# 
#   )  %>%
#   add_lines(x = ~E, y = ~phi_MB, type = 'scatter', mode = 'lines', showlegend=T, name=TeX("\\phi_{MB}")) %>%
#   config(mathjax = "cdn") %>%
#   add_bars(wap_f3, x = ~E, y = ~, type = 'scatter', mode = 'lines', showlegend=T, name=~shell, yaxis="y2")

# ### DO SUMS ETC FOR NEEC RATE


## Graph 73GE

# plot_ly(wap_f3, x = ~Ee, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% #, marker = list(color = label_colours)) %>%
#   layout(yaxis = list(type = "log", exponentformat = "E", title = "73GE NEEC Resonance Strength / beV"), #range = c(1E-10, max(wap_filtered$S))),
#          xaxis = list(tickangle = -45, title="Electron Impact Energy / keV", range = c(min(E), max(E))),
#          legend=list(name=~shell, showlegend=TRUE)
#          # yaxis2 = list(
#          #   tickfont = list(color = "red"),
#          #   overlaying = "y",
#          #   side = "right",
#          #   title = "MB" ,
#          #   showlegend = F)
#          
#   )  %>%
#   add_lines(x = ~E, y = ~phi_MB, type = 'scatter', mode = 'lines', showlegend=T, name="Phi_MB") %>%
#   add_lines(x = ~E, y = ~phi_FD, type = 'scatter', mode = 'lines', showlegend=T, name="Phi_FD")

## Graph
plot_ly(data=MBFD, x = ~E, y = ~phi_FD, type = 'scatter', mode="lines", showlegend=TRUE, name=TeX("\\phi_{FD}")) %>% #, marker = list(color = label_colours)) %>%
  layout(yaxis = list(range=c(log10(1E27), log10(1E34)), side="left",type='log', exponentformat = "E", title = TeX("\\phi(E) cm^{-2}s^{-1}eV^{-1}")), #range = c(1E-10, max(wap_filtered$S))),
         xaxis = list(showgrid=T, tickangle = -45, title="Electron Impact Energy  (keV)", range = c(min(E), max(E)), dtick=1),
         legend=list(showlegend=TRUE, x = 1.1, y = 0.9),
         yaxis2 = list(
           #tickfont = list(color = "red"),
           overlaying = "y",
           side = "right",
           title = "Resonance Strength (beV)" ,
           showlegend = F,
           type="log",
           exponentformat="E",
           showgrid=F
           
         )
         #range = c(0, Ebeam_MeV))
         #egend = list(x = 0.1, y = 0.9)
         
  )  %>%
  add_lines(x = ~E, y = ~phi_MB, type = 'scatter', mode = 'lines', showlegend=T, name=TeX("\\phi_{MB}")) %>%
  config(mathjax = "cdn") %>%
  add_bars(data=wap_f3, x = ~Ee, y = ~S, type = 'bar', mode = 'lines', showlegend=T, name=~shell, yaxis="y2")

#subplot(p1,p2, nrow=2, titleX=T, titleY=T)

################# Double Axis Plot #############################################
# plot_ly(Cand, x = ~x, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% #, marker = list(color = label_colours)) %>%
#   layout(plot_bgcolor='#e5ecf6',
#          yaxis = list(type = "log", exponentformat = "E", title = "NEEC Resonance Strength (beV)" , side="right"  ), #range = c(1E-10, max(wap_filtered$S))),
#          yaxis2 = list(
#            tickfont = list(color = "red"),
#            overlaying = "y",
#            side = "left",
#            title = "Energy of Ion Inside Target (MeV)" ,
#            showlegend = F,
#            range = c(0, Ebeam_MeV)),
#          xaxis = list(tickangle = -45, title="Target Depth (um)", range = c(0, 80), dtick = 10, showgrid=T),
#          legend=list(name=~shell, showlegend=TRUE)
#   ) %>%
#   add_lines(x = ~E_explore$x, y = ~E_explore$Eb, type = 'scatter', mode = 'lines', color=list("black"), yaxis="y2", showlegend=T, name="Ion Energy ") %>%
#   config(mathjax = "cdn") %>%
#   layout(legend= list(x=0.8, y=0.9))

################## LOOKING AT THE LIMITING DENSITIE ACCORDING TO Te ############

## Boltzmann constant
kb <- 8.617333262E-5 ##eVK^-1

## set T from Te ...eV
#Temp <- Te/kb
# Temp <- 4000 ##k
# Te <- Temp*kb
#Te <- 100000
Te <- 50
ne <- 1E24
Temp <- Te/kb

## deBroglie wavelength
lam_DB <- (1.0788E3)/(sqrt(Temp))  #### angstrom
lam_DB
##Set ne ... m^-3
#ne <- 1E27
ne_A <- ne*1E-30
##distance between electrons ...Ang
distance <- ne_A^(-(1/3))
distance

### Write a function that finds the limiting number density for each temperature
### SET: 
#Te <- 2500 #eV
ne_cm <- ne #cm^-3

ne_m <- ne_cm * 1E6 #...m^-3
### Look at degeneracy parameter
#Fermi energy
Ef <- (((hbar^2)/(2*me_kg))*((3*(pi^2)*ne_m))^(2/3))/e
#Ef <- ((hbar_eV^2)/(2*me_cms))*((3*(pi^2)*(ne*1E-6))^(2/3))
Ef

Theta <- Te/Ef
Theta



