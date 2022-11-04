                  library(tidyverse)
                  library(stringr)
                  library(plotly)
                  library(XML)
                  library(processx)
                  library(pbapply)
                  # install.packages("tidyverse")
                  # install.packages("stringr")
                  # install.packages("plotly")
                  # install.packages("XML")
                  # install.packages("processx")
                  # install.packages("pbapply")
                  
                  #### MAXIMUM AMOUNT OF ENERGY WE CAN PUT IN
                  Emaximum <- 4000 #keV
                 
                  #### FUNDAMENTAL CONSTANTS ####
                  e <- 1.602E-19 #C
                  mp <- 931.5E6 #eV/c^2
                  mn  <- 939.57E6 #eV/c^2
                  me <- 0.51099895000E6 #eV/c^2
                  me_kg <- 9.1093837E-31 #kg
                  hbar <- 1.054571817E-34	#J.s
                  hbar_eV <- 6.582119569E-16 #eVs
                  h <- hbar*2*pi
                  h_eV <- 4.135667696E-15 #eVs^-1
                  #A.U. to SI conversion
                  Eh <- 4.359744722E-18 ### one hartree is that many Joules
                  Eh_eV <- 27.211386 ## one hartree is that many electron volts
                  
                  #SI to a.u conversion
                  time_conv <- 2.418884326E-17 # this many seconds is an atomic unit ...hbar/Eh
                  length_conv <- 5.291772E-11 # this many m is one atomic unit ...length of one Bohr radius a0
                  barn <- 1E-28 #m^2
                  me.mp <- 1/(1836.15267343)
                  mp.me <- 1/me.mp
                  c <- 299792458 #ms^-1
                  me_cms <- me/((100*c)^2) #eV/cm-2s2
                  
                  ####### THIS ALL SHOULD ONLY BE RUN ONCE AT THE BEGINNING AFTER CLEARING WORKSPACE #######
                  source("tidy_functions2.R")
                  source("DecisionFunctions.R")
                  source("Microscopic/ICC_Eval.R")
                  ###### CALCULATION FUNCTIONS ######
                  source("Microscopic/WeisskopfCalculator.R")
                  source("Microscopic/PDB_ResonanceStrength_Calculator.R")
                  source("Macroscopic/Rates_EBIT.R")
                  source("Macroscopic/Rates_Plasma.R")
        
                
                  #### Read Data ####
                  ##read the gamma info##
                  ### Start by binding additions for B values that were added by hand ###
                  
                  #wap0 <- read_csv("App_Data/All_RatesCalcd_Viva.csv")

                  Gamma_Bvals_0 <- read_csv("Dependencies/ND_Gamma_Cont.csv")
                  Gamma_Bvals_Additions <- read_csv("Dependencies/Additional_Bvalues.csv")
                  Gamma_Bvals <- bind_rows(Gamma_Bvals_0, Gamma_Bvals_Additions)
                  Bvals_New <- read_csv("Dependencies/ND_Gamma_NEW.csv")

                  ##Electron binding Energies
                  Ebinds <- read.csv("Dependencies/Ebinds/NIST_Ebinds.csv", header=TRUE, fill=TRUE)
                  # 
                  ## NNDC  GS and isomers ... all good for "J"
                  #ND <- read_csv("Dependencies/ND.csv")
                
                 
                  ### FOR WHEN CHANGES TO tidy_Jpi are made or database reread ###
                  # ND <- tidy_Jpi(ND)
                  #ND <- unique(ND)
                  #ND <- filter(ND, Multiple_Jpi == FALSE) 
                  # write.table(ND, file="Dependencies/ND_TidyJpi.csv", append=T, row.names=F, col.names=T,  sep=",")
                  ## Read the tidy nuclear data
                  ND <- read_csv("Dependencies/ND_TidyJpi.csv")
                  
                  # 
                  # ### GROUND STATES ###
                  GS <- filter(ND, E==0, Multiple_Jpi == FALSE, !is.na(Thalf))
                  #GS <- read_csv("GS.csv")
                  GS2 <- rename(GS, Ji = J, Thi = Thalf)
                  for(i in 1:length(GS2$Thi)){
                    print(i)
                    #i <- 3
                    GS2$Thi_numeric[i] <- as_Numeric_Th(GS2$Thi[i])
                    #as_Numeric_Th(GS2$Thi[i])
                  }
                  # ####### DONT WRITE GS2 FOR NOW JUST LEAVE GS2 AS IS AS IT WORKS AND NOT TOO SLOW UP TO HERE #######
                  # 
                  
                  ### Alternative List of isomers
                  # ND$Thi_numeric <- pblapply(ND$Thalf, as_Numeric_Th)
                  # NDi <- filter(ND, Thi_numeric > 1E-9, E != 0)
                  # 
                  # NDi2 <- distinct(NDi, AX, E, J, .keep_all=TRUE)
                  
                  
                  
                  ### Isomer ###
                  #isomer <- read_csv("Dependencies/isomer.csv")
                  source("isoTidy.R")
                  # #isomer <- read_csv("Dependencies/isomer.csv")
                  isomer <- tidy_Jpi(isomer) %>% filter(Multiple_Jpi == FALSE) #%>% rename(Ji = J)
                  for(i in 1:nrow(isomer)){
                      if(  !str_detect(isomer$J[i], "[\\+\\-]")  ){
                          isomer$J[i] = paste0(isomer$J[i], "-")
                      }
                  }


                  # unlink("Dependencies/isomer.csv")
                  # write.table(isomer, file="Dependencies/isomer.csv", append=T, row.names=F, col.names=T,  sep=",")
                  #write.table(ND, file="ND_TidyJpi.csv", append=T, row.names=F, col.names=T,  sep=",")
                  
                  ## Handy to have iso_AX for the app
                  #wap_Ni <- read_csv("Outputs/Wins_all_plasma_Ni.csv", col_types = cols( Ji = col_character(), ICC = col_double(), ICC_tot = col_double(), S = col_double(), S_tot = col_double(), N_neec_pw_Stot = col_double(), N_neec_pw = col_double()   ))
                  #iso_AX <- unique(filter(wap_Ni, !is.na(S), S > 0)$AX)
                  
                  
                  # #### SOME NEW STUFF FOR ISOMER COMPARISON
                  # ND$Thalf_numeric <- pblapply(ND$Thalf, as_Numeric_Th)
                  # ND_iso <- filter(ND, E>0, )
                  
                  # Mo <- read_csv("App/All_RatesCalcd_Viva.csv")
                  # Mo <- filter(Mo, AX=="93MO") %>%
                  #     arrange(Ef)
                  # 
                  # write.table(Mo, file="Mo_Addition.csv", append=T, row.names=F, col.names=T,  sep=",")
                  # 
                  ## STATE number of electrons on ion (occ) (after capture) 
                  #n <- 1
#source("Macroscopic/All_Rates_Compute.R")
                  
                 
                 
                  # ## CHUNK ND
                  # # 
                  # #                   ## Then regroup
                  # #                   #for(i in 1:length(group))
                  #                   args <- commandArgs(TRUE)
                  #                   m1 <- as.numeric(args[1])
                  #                   m2 <- as.numeric(args[2])
                  #                   source("Macroscopic/All_Rates_Compute.R")
                  # #                   # m1 <- 1
                  # #                   # m2 <- 298
                  # #                    # ND <- ND %>% select(-DE, -DThalf, -Itf)
                  # #                   # ND <- ND %>% filter(M >= m1 & M <= m2)
                  # 
                  # 

    # #            # for(n in 79:89){
    # #              n <- 1
    # #            # #      #n <- 47
    # #            # # #
    #              args <- commandArgs(TRUE)
    #              n <- as.numeric(args[1])
    #              # n<-100
    # #            # #
    # #            #
    #              Name_AnN0_Energetic <- sprintf("Outputs/MASTER/A%dN0_Plasma_Emax.csv", n)
    #              Name_AnNi_Energetic <- sprintf("Outputs/MASTER/A%dNi_Plasma_Emax.csv", n)
    # 
    #              Name_AnN0_Rates <- sprintf("Outputs/MASTER/A%dN0_Plasma_Rates.csv", n)
    #              Name_AnNi_Rates <- sprintf("Outputs/MASTER/A%dNi_Plasma_Rates.csv", n)
    # #            # #
    # #            # #
    # #            # #     # Name_AnN0_Energetic <- sprintf("Outputs/2021/Plasma/Energetics/A%dN0_Plasma_Emax.csv", n)
    # #            # #     # Name_AnNi_Energetic <- sprintf("Outputs/2021/Plasma/Energetics/A%dNi_Plasma_Emax.csv", n)
    # #            # #     #
    # #            # #     # Name_AnN0_Rates <- sprintf("Outputs/2021/Plasma/Rates/A%dN0_Plasma_Rates.csv", n)
    # #            # #     # Name_AnNi_Rates <- sprintf("Outputs/2021/Plasma/Rates/A%dNi_Plasma_Rates.csv", n)
    # #            # #     #
    # #            # #      #### READING EBIT Energetic results applyig plasma analysis
    # #            # #      # Name_AnN0 <- sprintf("Outputs/2021/R-EBIT_30keV/Energetics/A%dN0_R-EBIT_30keV.csv", n)
    # #            # #      # Name_AnNi <- sprintf("Outputs/2021/R-EBIT_30keV/Energetics/A%dNi_R-EBIT_30keV.csv", n)
    # #            # #
    # #            # #
    # #            # #      # Name_AnN0_Rates <- sprintf("Outputs/2021/Plasma/Rates/A%dN0_Plasma_Rates.csv", n)
    # #            # #      # Name_AnNi_Rates <- sprintf("Outputs/2021/Plasma/Rates/A%dNi_Plasma_Rates.csv", n)
    # #            # #    ### is using Vi then i is zero indexed to imply the occupation before capture
    # #            # #
    # #            # #    ### FILENAME PROTOTYPE: AnNi_*Approach*_*Facility*    ... for n electron ion, isomeric ###
    # #            # #    #Name_AnN0 <- sprintf("Outputs/2020/A%dN0_LT15MeV_EBITRates.csv", n)
    # #            # #    #Name_AnNi <- sprintf("Outputs/2020/A%dNi_LT15MeV_EBITRates.csv", n)
    # #            # #
    # #            # #    ## Allocate key column
    # #            # #    #ND <- ND %>% mutate(Key=c(1:length(ND$E)), I1tf=logical(1), I2tf=logical(1), I3tf=logical(1), I4tf=logical(1), I5tf=logical(1), I6tf=logical(1), I7tf=logical(1), I8tf=logical(1), I9tf=logical(1), I10tf=logical(1), I11tf=logical(1))
    # #            # #    # WbeamI_k_tf <- tibble(Key=c(1:length(ND$E)), I1tf=logical(1), I2tf=logical(1), I3tf=logical(1), I4tf=logical(1), I5tf=logical(1), I6tf=logical(1), I7tf=logical(1), I8tf=logical(1), I9tf=logical(1), I10tf=logical(1), I11tf=logical(1))
    # #            # # #
    #               ## Should we do energy loops
    #               Energetic <- F
    #               ## Should we do the wrangling after
    #               Wrangling <- T
    # 
    #              # source("AnN0_MapReplace.R")
    #               source("AnNi_MapReplace.R")
    # #            # # #
    # #            # # #
               # #     source("AnN0_MapReplace.R")
               # # #    #source("AnN0.R")
               #
               #    Energetic <- F
               #    ## Should we do the wrangling after
               #    Wrangling <- T
               #     source("AnNi_MapReplace.R")
    #            #    #source("AnNi.R")
    # #}
               #    #source("EexcSearch.R")

                ## Testing
                # n <- 1
                # level <- ND[414836,]
                # eb <- filter(Ebinds, Z==level$Z) %>% filter(Occ==n)
                # eb$IE

                ################# HAVING A QUICK LOOK AT B VALUE TRENDS IN Gamma_Bvals


                
                
                
                
                
                
                
                ## Write a certain nuclide 
                # Mo <- ND %>% filter(AX == '93MO')
                # write.table(Mo, file="~/Documents/current/Candidate/Specific_Candidates/93Mo.csv", append=T, row.names=F, col.names=T,  sep=",")
                
                
                #### Goodness of experimental conditions columns ####
                
                # #read csv and mutate accordingly
                # isolde <- read_csv("Dependencies/isolde_Beams.csv")
                # caribu <- read_csv("Dependencies/caribu_Beams.csv")
                # ND_Beams <- read_csv("Outputs/AnNi/A0Ni_beam_001_298")
                # 
                # isolde$A <- as.numeric(str_extract(isolde$Anum, "[0-9]"))
                # isolde$GorI <- str_extract(isolde$Anum, "[aA-zZ]+")
                # 
                # caribu$A <- as.numeric(str_extract(caribu$Nuclide, "[0-9]"))
                # caribu$El <- str_extract(caribu$Nuclide, "[aA-zZ]+")
                # 
                # #add facility tf column
                # isomer %>% mutate(isolde_tf=logical(1), caribu_tf=logical(1))
                # 
                # 
                # ## Do some plots
                # A0Ni <- read_csv("Outputs/AnNi/A0Ni_beam_001_298.csv")
                # A1Ni <- read_csv("Outputs/AnNi/A1Ni_beam_001_298.csv")
                # A2Ni <- read_csv("Outputs/AnNi/A2Ni_beam_001_298.csv")
                # A3Ni <- read_csv("Outputs/AnNi/A3Ni_beam_001_298.csv")
                # A4Ni <- read_csv("Outputs/AnNi/A4Ni_beam_001_298.csv")
                # A5Ni <- read_csv("Outputs/AnNi/A5Ni_beam_001_298.csv")
                # A6Ni <- read_csv("Outputs/AnNi/A6Ni_beam_001_298.csv")
                # A7Ni <- read_csv("Outputs/AnNi/A7Ni_beam_001_298.csv")
                # A8Ni <- read_csv("Outputs/AnNi/A8Ni_beam_001_298.csv")
                # A9Ni <- read_csv("Outputs/AnNi/A9Ni_beam_001_298.csv")
                # 
                # par(mfrow=c(2, 3))
                # 
                # plot((A0Ni$M-A0Ni$Z), A0Ni$Z, col='red', pch=".", xlab="N", ylab="Z")
                # title("n=0")
                # plot((A1Ni$M-A1Ni$Z), A1Ni$Z, col='orange', pch=".", xlab="N", ylab="Z")
                # title("n=1")
                # plot((A2Ni$M-A2Ni$Z), A2Ni$Z, col='yellow', pch=".", xlab="N", ylab="Z")
                # title("n=2")
                # plot((A3Ni$M-A3Ni$Z), A3Ni$Z, col='green', pch=".", xlab="N", ylab="Z")
                # title("n=3")
                # plot((A4Ni$M-A4Ni$Z), A4Ni$Z, col='blue', pch=".", xlab="N", ylab="Z")
                # title("n=4")
                # plot((A5Ni$M-A5Ni$Z), A5Ni$Z, col='black', pch=".", xlab="N", ylab="Z")
                # title("n=5")
                # points((A6Ni$M-A6Ni$Z), A6Ni$Z, col='violet', pch=".", xlab="N", ylab="Z")
                # points((A7Ni$M-A7Ni$Z), A7Ni$Z, col='orange', pch=".", xlab="N", ylab="Z")
                # points((A8Ni$M-A8Ni$Z), A8Ni$Z, col='orange', pch=".", xlab="N", ylab="Z")
                # 
                # 
                # 
                # 
                
                
                
                
                
                
                
                
                 
