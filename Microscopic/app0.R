
library(rsconnect)
library(shiny)
library(tidyverse)
library(plotly)
library(DT)  ## for data tables



#### MAIN DATA TABLES
#ND <- read_csv("OUTPUTS/ND_Gamma.csv")
library(rsconnect)
library(shiny)
library(tidyverse)
library(plotly)

## set working directory as Microscopic/Scaling directory
setwd("../")

source("Wins.R")

## Read required data tables ... these are the newest tables for conversion of quantuym labels
n_index <- read_csv("Microscopic/n_to_ICshell_Conversion.csv")
Ebinds <- read.csv("Dependencies/Ebinds/NIST_Ebinds.csv", header=TRUE, fill=TRUE)
#LS_Jat_Conv <- read.csv("Microscopic/elec_Jconv.csv")
## Add a dirac angular momentum quantum number variable
#ICC_multipliers <- read.csv("Microscopic/ICC_Multipliers.csv")



# ##### Wrangle ####
# dat <- read.table("Microscopic/Scaling/mo93.neec-data-NOHEAD.3wsu.dat", sep="")
# Abinit <- as_tibble(dat)
# colnames(Abinit) <- c("q", "n", "kappa", "J_final", "Eres_eV", "Y_neec_SI", "S_neec_beV")

#### Read database Mo ######
Mo <- read_csv("Mo_Additions_Rates2.csv")
Q_eV <- 4.85E3 #eV
Q <- 4.85

Abinit <- read_csv("Microscopic/Abinit_Apr2022.csv")

########## RUN FROM HERE WHEN DOING RATES COMPARISON GRAPHS
###### Adjust Mo PDB Resonances for non gamma scaled   ######
Mo <- read_csv("Mo_Additions_Rates2.csv")
Mo <- filter(Mo, !is.na(S))
Mo <- mutate(Mo, fq=double(1))
for(i in 1:length(Mo$S)){
  Mo$S[i] <- Mo$S[i] * 1/(Mo$Gamma_Scale_Factor[i])
}
#Check against Abinit GS
Abinit_GS <- filter(Abinit, GS==TRUE)


N_out <- read_csv("/home/ben/physnp/Candidate/neec-candidate/Optims/N_out_93mMo_4.85keV_AllFLYCHK.csv")
#N_out <- read_csv("Optims/N_out_57Fe_14keV_AllFLYCHK.csv")

N_out$Te <- as.double(N_out$Te) / 1000
N_out$ne <- as.double(N_out$ne)
N_out$N_count <- as.double(N_out$N_count)
N_out <- filter(N_out, N_count > 1E-20)


### See if the optimum temperature deviates accross density

##Create array of optimum temps
Temp_Optims <- tibble(Temp = double(length(FLYCHK_ne)), ne=FLYCHK_ne, Rate=double(length(FLYCHK_ne)) )

## factor temperature sets for each density
for(i in 1:length(FLYCHK_ne)){
  N_out2 <- filter(N_out, ne == FLYCHK_ne[i] )
  Temp_Optims$Temp[i] <- N_out2$Te[which.max(N_out2$N_count)]
  Temp_Optims$Rate[i] <- max(N_out2$N_count)
}







####################################################################################################
### UI #############################################################################################
####################################################################################################
ui <- fluidPage(
  ## Select Variable Inputs
  selectInput("elec_d", "Select Electron Number Density ", choices=FLYCHK_ne, selected="1E20"),
  plotlyOutput("plot1"),
  plotlyOutput("plot2")
  
)

####################################################################################################
## SERVER ##########################################################################################
####################################################################################################
server <- function(input, output) {
  
  Temp_reactive <- reactive({
    ne_set <- as.double(input$elec_d)
    ##### RUN over the NEEC spectrum energies and decide on best Te for Ee = Te #######
    Temps <- tibble(Te=(FLYCHK_Te*1E-3), Rtot=double(1), Rtot_psum = double(1), Rtot_PDBsum=double(1))
    for(j in 1:length(Temps$Te)){
      #### At this point S_sum is an order of magnitude too large
      #### Using Stot from alpha_tot = 4.63E5
      ######### Get average charge state #############
      ne <- ne_set ### e cm^-3
      #Te <- 1.441
      #<-22
      Te <- Temps$Te[j] ###keV
      Te_eV <- Te * 1E3
      Z_name <- sprintf("/home/ben/physnp/Candidate/neec-candidate/Dependencies/CSD/CSD_Ave_Data/Z%003.f_AVE_CS.csv", 42)
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
      CSD_File <- sprintf("/home/ben/physnp/Candidate/neec-candidate/Dependencies/CSD/CSD_fq_Data/CSD_%d.csv", 42)
      CSD <- read_csv(CSD_File)
      CSD_neTe <- filter(CSD, ne_cm3 == ne_choice, Te_eV == Te_choice)
      
      #### We read the charge state of the CS chosen by the closest Te and ne selected by the laser params
      Z_ave <- filter(Z_dat, Te == Te_choice, ne == ne_choice)$CS
      #fq <- CSD_neTe$frac[which(CSD_neTe$q==(Wins_all$CS+1))]
      #if(Wins_all$Z >)
      
      
      ###### Compare with plasma yields #####
      ## Literature Yield .... Te = 2 keV, ne = 1E23, Rneec = 5 s^-1
      ## S_tot Yield
      # ne <- 1E23 ### e cm^-3
      # Te <- 2 ###keV
      alpha_tot <- 4.63E5
      #compute_S <- function(Jd, Jat,Ji,IC, Ar, Etrans_keV, Ebind_keV){
      
      ## Volume and number of ions
      Rp<- 40E-4 #cm
      Vp <- (4/3) * pi * (Rp^3) #cm^3
      ni <- (ne / Z_ave) 
      f_iso <- 1E-5
      Ni <- ni * Vp * f_iso
      
      
      ### assign fq based on current {Te,ne} ###
      for(i in 1:length(Abinit$q)){ 
        Abinit$fq[i] <- CSD_neTe$frac[which(CSD_neTe$q==(Abinit$q[i]))]
      }
      for(i in 1:length(Mo$CS)){ 
        Mo$fq[i] <- CSD_neTe$frac[which(CSD_neTe$q==(Mo$CS[i]+1))]
      }
      
      
      ### Calculate neec rate for each capture channel
      for(i in 1:length(Abinit$q)){
        Abinit$R_neec_p[i] <- ne * MB(Abinit$Eres_eV[i], 1E3*Te) *  sqrt(1 - (1/((1+(Abinit$Eres_eV[i]) / me))^2))*100*c   *  Abinit$S_neec_beV[i] * 1E-24 * (2/(pi)) *  Abinit$fq[i] 
      }  
      for(i in 1:length(Mo$S)){
        Mo$Rneec_Pl[i] <- ne * MB(Mo$Ee[i]*1E3, 1E3*Te) *  sqrt(1 - (1/((1+(Mo$Ee[i]*1E3) / me))^2))*100*c   *  Mo$S[i] * 1E-24 * (2/(pi)) *  Mo$fq[i] 
      }
      R_neec_pTot <- sum(Abinit$R_neec_p)
      R_neec_PDBTot <-sum(Mo$Rneec_Pl)
      
      
      #### What is the best energy that gives the best sum rate
      # set parameters for sum rate 
      #ne <- 1E24 ### e cm^-3
      #Te <- 1.441 ###keV
      # Alpha tot
      #Vi_chose <- mean(Abinit$Vi_keV)
      Vi_chose <-  mean(filter(ICC_Multipliers, Z==42)$Ebind)
      Q<-4.85
      alpha_tot <- 4.63E5
      Ee_chose <- Q-Vi_chose
      #E_max <- Q - mean(filter(ICC_Multipliers, Z==42)$Ebind)
      #E_max <- Q - max(Abinit$Vi_keV)
      S_beV <- compute_S(8.5, mean(Mo$Jat), 10.5, alpha_tot, Mo$Ar[10], Q, Vi_chose)
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
      
      R_neec_tot <- S * (2/(pi)) * F_ETe * N_res * ne * v_max 
      
      Temps$Rtot[j] <- R_neec_tot
      Temps$Rtot_psum[j] <- R_neec_pTot
      Temps$Rtot_PDBsum[j] <- R_neec_PDBTot
      
      # R_neec_tot
      # R_neec_pTot
      
    }
    return(Temps)
  })

  output$plot1 <- renderPlotly({
    
    Temps_plot <- Temp_reactive()
    
    plot_ly(Temps_plot, x = ~Te,y= ~Rtot_PDBsum, type='scatter', mode="markers+lines", name="Rtot_PDBsum", showlegend=TRUE) %>% #, marker = list(color = label_colours)) %>%
      layout(yaxis = list(  title = "NEEC Rate Per Ion (s^(-1))"), #range = c(1E-10, max(wap_filtered$S))),
             xaxis = list(tickangle = -45, title="Temp [keV]", range=c(0.1,30)) #%>%
             #legend=list(name=~shell, showlegend=TRUE) %>%

             #   yaxis2 = list(
             #   tickfont = list(color = "red"),
             #   overlaying = "y",
             #   side = "right",
             #   title = "MB" ,
             #   showlegend = F)
      )  %>%
      add_trace(x= ~Te,  y = ~Rtot, type = 'scatter', mode='markers+lines', name="R_Alpha_tot") %>%
      add_trace(x= ~Te, y= ~Rtot_psum, type='scatter', mode="markers+lines", name = "Rtot_PalffyWu")
    
    
    # plot_ly(Abinit, x = ~Ee, y = ~S, type = 'bar', showlegend=TRUE, name="Resonance Strength (beV)") %>% #, marker = list(color = label_colours)) %>%
    #   layout(yaxis = list(  title = "NEEC Rate Per Ion (s^(-1))"), #range = c(1E-10, max(wap_filtered$S))),
    #          xaxis = list(tickangle = -45, title="Temp [keV]", range=c(0.1,30)) #%>%
    #          #legend=list(name=~shell, showlegend=TRUE) %>%
    #          
    #          #   yaxis2 = list(
    #          #   tickfont = list(color = "red"),
    #          #   overlaying = "y",
    #          #   side = "right",
    #          #   title = "MB" ,
    #          #   showlegend = F)
    #   )  %>%
    #   add_lines( x= ~Temps_plot$Te, y= ~Temps_plot$Rtot_psum, type='scatter', mode="markers+lines", name="Rtot_psum") %>%
    #   add_lines( x= ~Temps_plot$Te, y= ~Temps_plot$Rtot_PDBsum, type='scatter', mode="markers+lines", name = "Rtot_PDBsum") 
    
    #   
      
    
  })  
  
  output$plot2 <- renderPlotly({
  plot_ly(N_out, x = ~Te, y= ~ne,  z = ~N_count, type="scatter3d", mode = 'markers', marker=list(size=3)) %>%
    layout(scene=list(
      xaxis=list(type = "log", title="Temperature (keV)", exponentformat = "E"),
      yaxis=list(type = "log", title="Electron number density (cm^(-3))", exponentformat = "E"),
      zaxis=list(type = "log", title="Total NEEC Rate (s^-1)", exponentformat = "E")
    )
    )
  
   })  
  
  
}

shinyApp(ui = ui, server = server)
