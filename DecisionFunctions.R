#### WIN DECISION FUNCTIONS ####


### CONSTRAINTS ###
#Beam must be less than 15MeV/nucleon
#beamCons <- 15
#Effective electron kinetic energy constraint
#Ee <- beamCons*(1/1836)*1000
### CURRENTLY SET AT BEGINNING OF Wins.R script
Ee <- Emaximum #keV

#### Modified winbeamtf function to allow for New atomic data structure
winBeamtf <- function(level, Einit, n) {
  # check not NA
  if(is.na(level$E) || is.na(Einit) || (level$Z>103) ) return(NA)
  if(  (level$El=="RF") && (n<90) && (n>78) ) return(NA)
  Q <- level$E - Einit
  ### **  # Ebinds column to read
  ### Ebinds row to read
  ## Select the occupancy to read based on fed in value of n
  Ebind_level <- filter(Ebinds, Z==level$Z) %>% filter(Occ==n)
  if(n>level$Z) return(NA)
  Vi <- Ebind_level$IE
  if(is.na(Vi)) return(NA)
  Emax <- Vi + Ee
  if(  (Vi <= Q) && (Q <= Emax) ) { return(TRUE) }
  else return(FALSE)
}

#### Modified winbeamtf function to allow for New atomic data structure
######Further modified to make inputs single
winBeamtf_SV <- function(levelZ, levelE, Einit, n) {
  # check not NA
  if(is.na(levelE) || is.na(Einit) || (levelZ>103) ) return(NA)
  Q <- levelE - Einit
  ### **  # Ebinds column to read
  ### Ebinds row to read
  ## Select the occupancy to read based on fed in value of n
  Ebind_level <- filter(Ebinds, Z==levelZ) %>% filter(Occ==n)
  if(n>levelZ) return(NA)
  Vi <- Ebind_level$IE
  if(is.na(Vi)) return(NA)
  Emax <- Vi + Ee
  if(  (Vi <= Q) && (Q <= Emax) ) { return(TRUE) }
  else return(FALSE)
}


## Vec_ ...  functions that take single value inputs and are compatible with pmap
Vec_winExctf <- function(E, Einit, Qmin, Qmax) {
  if(is.na(E) || is.na(Einit)) return(NA)
  ## Evalutate nuclear excitation energy based on input
  Q <- E - Einit
  if(  (Qmin <= Q) && (Q <= Qmax)  ) { return(TRUE) }
  else return(FALSE)
}

#pmap(list(1,2,3,4), Vec_winExctf)
####WORKS RETURNS FALSE###



# ## Win beam TRUE/FALSE based on Ee constraints, n electron spit-orbit shell (no splitting)
# winBeamtf <- function(level, Einit, n) {
#   # check not NA
#   if(is.na(level$E) || is.na(Einit)) return(NA)
#   Q <- level$E - Einit
#   # Ebinds column to read
#   ecol <- n+2
#   Vi <- Ebinds[level$Z, ecol]
#   if(is.na(Vi)) return(NA)
#   Emax <- Vi + Ee
#   if(  (Vi <= Q) && (Q <= Emax) ) { return(TRUE) }
#   else return(FALSE)
# }

## Beam TF return isomer n-tuple of wins for that nuclide level
winBeamtfIso <- function(level, Einit, ind, n) {
  if(is.na(level$E) || is.na(Einit) || (levelZ>103)) return(NA)
  Q <- level$E - Einit
  #old ... ecol <- n+2
  #old ... Vi <- Ebinds[level$Z, ecol]
  
  ### **  # Ebinds column to read
  ### Ebinds row to read
  ## Select the occupancy to read based on fed in value of n
  Ebind_level <- filter(Ebinds, Z==level$Z) %>% filter(Occ==n)
  if(n>level$Z) return(NA)
  Vi <- Ebind_level$IE
  
  if(is.na(Vi)) return(NA)
  #Ei <- Ebinds$K[level$Z]
  Emax <- Vi + Ee
  #if(  (Ei <= Et) && (Et <= Emax) ) { return(ind) }
  if(  (Vi <= Q) && (Q <= Emax) ) { return(isomer$Key[ind]) }
  else return("-1")
}


############# LABEL FIRST OBSERVATIONS ###############
# ##TODO .... allocate a labelling to number which Ef
# ## How to allocate the FIRST OBSERVATION
# for(i in 1:length(Wap$Ef)){
#   
#   if(  (i+1 > length(Wap$Ef))  ){
#     if( (Wap$first_Ef[i-1] == TRUE)  ){
#       Wap$first_Ef[i] <- TRUE
#       next
#     } else {
#       Wap$first_Ef[i] <- FALSE
#       next
#     }
#   }
#   
#   #if it is the same as the last one and has the same Jf, Type
#   if( (switch == 1) && (Wap$Ef[i] == Wap$Ef[i+1]) && (Wap$Jf[i] == Wap$Jf[i+1]) && (Wap$Type[i] == Wap$Type[i+1])  ) {
#     Wap$first_Ef[i+1] <- TRUE
#   }
#   
#   ##If its not the same as the last one and has the same Jf, Type
#   if( (switch == 1) &&  (Wap$Ef[i] != Wap$Ef[i+1]) && (Wap$Jf[i] == Wap$Jf[i+1]) && (Wap$Type[i] == Wap$Type[i+1])   ) {
#     Wap$first_Ef[i+1] <- FALSE
#     #store the index of the last first obs
#     last_Ef <- i
#     switch <- 2
#     
#     # ## All i's are false now untill a new Ef comes along
#     # if((Wap$Ef[i] != Wap$Ef[i+1])){
#     #   ##RESTART
#     #
#     # }
#     #
#   }
#   
#   # if( (switch == 1) && (Wap$Jf[i] != Wap$Jf[i+1]) && (Wap$Type[i] != Wap$Type[i+1])  ) {
#   #   switch <- 1
#   # }
#   
#   if( (switch == 2) &&  (Wap$Ef[i] == Wap$Ef[i+1]) && (i != length(Wap$Ef)) ){
#     Wap$first_Ef[i] <- FALSE
#     Wap$first_Ef[i+1] <- FALSE
#     switch <- 2
#   }
#   if( (switch == 2) && (Wap$Ef[i] != Wap$Ef[i+1]) &&  (Wap$Jf[i] == Wap$Jf[i+1]) && (Wap$Type[i] == Wap$Type[i+1])  ) {
#     Wap$first_Ef[i+1] <- FALSE
#     switch <- 2
#   }
#   
#   if( (switch == 1) &&  (Wap$Ef[i] == Wap$Ef[i+1]) && (i != length(Wap$Ef)) ){
#     switch <- 1
#   }
#   if( (switch == 1) &&  (Wap$Ef[i] != Wap$Ef[i+1]) && (i != length(Wap$Ef)) ){
#     switch <- 2
#   }
#   
#   if( (switch == 2) && (Wap$Jf[i] != Wap$Jf[i+1]) && (Wap$Type[i] != Wap$Type[i+1])  ) {
#     Wap$first_Ef[i+1] <- TRUE
#     switch <- 1
#   }
#   
#   #within 1 keV?
#   #if(  (Wap$Ef[i] == (Wap$Ef[i] + 1)) && (Wap$Jf[i] )
#   
# }
# 



### Decision Based on Multipolarity ###
## INPUT an initial and final Jpi and return TYPE (M1 etc...)
# function prototype: MultiPOL("M1", "3/2", "1/2")
MultiPOL <- function(type, Jpi_i, Jpi_f) {
  
  if(is.null(Jpi_i) || is.null(Jpi_f)) (return(NA))
  ## condition based on M1 or E1 or E2 Declaration etc.
  #Declaration:
  condition <- type
  deltaJ <- mP$deltaJ[which(mP$Type==condition)]
  deltaPI <- mP$deltaPI[which(mP$Type==condition)]
  
  #return NA if NA entries
  if(  is.na(Jpi_i) || is.na(Jpi_f)  ){
    return(NA)
  }
  
  #if the first letter is an M then read off table for deltaPI
  
  #Extract the J and evaluate the actual values of the "eg. 9/2" strings
  Ji <- paste(unlist(str_extract_all(Jpi_i, "[0-9/]")), collapse = "")
  Ji <- eval(parse(text=Ji))  
  Jf <- paste(unlist(str_extract_all(Jpi_f, "[0-9/]")), collapse = "")
  Jf <- eval(parse(text=Jf))
  
  if(is.null(Ji) || is.null(Jf)) (return(NA))
  
  #Sort out the entries that are JUST + or -
  if(Jpi_i == "+" || Jpi_f == "+" || Jpi_i == "-" || Jpi_f == "-" ) {
    return(NA)
  }
  
  #Extract the parity of the level
  #if its "" then make sure its a "+"
  if(!str_detect(Jpi_i, "[+-]")) {pi_i <- "+"}
  pi_i <- paste(unlist(str_extract_all(Jpi_i, "[+-]")), collapse = "")
  
  if(!str_detect(Jpi_f, "[+-]")) {pi_f <- "+"}
  pi_f <- paste(unlist(str_extract_all(Jpi_i, "[+-]")), collapse = "")
  
  # check if the initial and final levels match the condition requirement
  if(   (abs(Jf - Ji) <= deltaJ)  && ((deltaPI == TRUE) && (pi_i != pi_f)) ||  ( (abs(Jf - Ji) == deltaJ)  &&  ( (deltaPI == FALSE) && (pi_i == pi_f) )  )  )     {
    return(TRUE)
  }
  
  else return(FALSE)
}

### Decision Based on Multipolarity ###
## INPUT an initial and final Jpi and return TYPE (M1 etc...)
# function prototype: MultiPOL("M1", "3/2", "1/2")

####### Labelling Functions Below ############
##############################################


mP <- read_csv("Dependencies/MultiPolarity_Selection.csv")

## Should input Ji and Jf and output a multipolarity type
MultiPOL2 <- function(Jpi_i, Jpi_f) {
  #print()
  if(is.null(Jpi_i) || is.null(Jpi_f)) (return(NA))
  ## condition based on M1 or E1 or E2 Declaration etc.
  #Declaration:
  #condition <- type
  #deltaJ <- mP$deltaJ[which(mP$Type==condition)]
  #deltaPI <- mP$deltaPI[which(mP$Type==condition)]
  
  # ## sort out the entry "3/2 5/2"
  # if (str_detect("3/2 5/2")) {
  #   return
  # }
  
  #return NA if NA entries
  if(  is.na(Jpi_i) || is.na(Jpi_f)  ){
    return(NA)
  }
  
  
  if(Jpi_i == Jpi_f) return("E0")
  
  
  #Sort out the entries that are JUST + or -
  if(Jpi_i == "+" || Jpi_f == "+" || Jpi_i == "-" || Jpi_f == "-" ) {
    return(NA)
  }
  
 
  
  ############ Handy for unstrnging strings paste(unlist( ... ), collapse = "")
  #Extract the J and evaluate the actual values of the "eg. 9/2" strings
  Ji <- paste(unlist(str_extract_all(Jpi_i, "[0-9/]")), collapse = "")
  Ji <- eval(parse(text=Ji))  
  Jf <- paste(unlist(str_extract_all(Jpi_f, "[0-9/]")), collapse = "")
  Jf <- eval(parse(text=Jf))
  
  if(is.null(Ji) || is.null(Jf)) (return(NA))
  
  if(Ji == Jf) return("E0")
  

  #Extract the parity of the level
  #if its "" then make sure its a "+"
  #if(!str_detect(Jpi_i, "\[+-\]")) {pi_i <- "+"}
  pi_i <- paste(unlist(str_extract_all(Jpi_i, "[+-]")), collapse = "")
  
  #if(!str_detect(Jpi_f, "[+-]")) {pi_f <- "+"}
  pi_f <- paste(unlist(str_extract_all(Jpi_f, "[+-]")), collapse = "")
  
  if(pi_i == pi_f) dPI <- FALSE
  if(pi_i != pi_f) dPI <- TRUE
  
  #make sure mP is read
  
  
  # # check if the initial and final levels match the condition requirement
  # if(   (abs(Jf - Ji) <= deltaJ)  && ((deltaPI == TRUE) && (pi_i != pi_f)) ||  ( (abs(Jf - Ji) == deltaJ)  &&  ( (deltaPI == FALSE) && (pi_i == pi_f) )  )  )     {
  #   return(TRUE)
  # }
  # 
  # else return(FALSE)
  ## Ensure the photon is an integer boson
  if( round(abs(Jf - Ji)) != abs(Jf - Ji)  ) return(NA)
  if(abs(Jf - Ji) > 5) return("Large_L")
  for(z in 1:12){
    # check if the initial and final levels match the condition requirement
    # if(  ( (abs(Jf - Ji) == mP$deltaJ[z])  && ((mP$deltaPI[z] == ) && (pi_i != pi_f))) ||  ( (abs(Jf - Ji) == mP$deltaJ[z])  &&  ( (mP$deltaPI[z] == FALSE) && (pi_i == pi_f) )  )  )     {
    if(   (abs(Jf - Ji) == mP$deltaJ[z]) && (dPI == mP$deltaPI[z])  ){
      return(mP$Type[z])
    }
  }
}

##### MAKE THi and THf numeric #####
## Need to sort out the ~symbol at the beginning of an input
as_Numeric_Th_iso <- function(Thi) {
    if(is.na(Thi)) {return(NA)}
    ## account for the stable isomer
    if(  str_detect(Thi, "STABLE")  ){return(-1)}
    ## discount funny symbol cases
    if(  str_detect(Thi, "@")  ){return(NA)}
    #else {return(0)}
    ## remove error and brackets
    if( str_detect(Thi, "[<>~\\*\\+\\-]*") ) {
      Thi <- str_remove_all(Thi, "[<>~\\*\\+\\-]*")
    }  
      
    if(  str_detect(Thi, "(\\(\\d+\\))|[<>~\\*\\+]")  ){
      #if(str_detect(Thi, "<>~*")) {
        Th <- str_remove_all(Thi, "(\\(\\d+\\))|[<>~\\*\\+]") 
        if( str_detect(Th, "[<>~\\*\\+]")  ){
          #if(str_detect(Thi, "]<>~*]")) {
          Th <- str_remove_all(Th, "[<>~\\*\\+]") 
          Th <- trimws(Th)
        }  
      #}
      }
    
    else if( str_detect(Thi, "[<>~\\*\\+]")  ){
      #if(str_detect(Thi, "<>~*")) {
      Th <- str_remove_all(Thi, "[<>~\\*\\+]") 
      Th <- trimws(Th)
    }  
    else (Th <- Thi)
   
    ## convert "  ns" to E-9
    if(  str_detect(Th, "(\\sns)")  ) {
      Th <- str_remove_all(Th, "(\\sns)")
      #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
      #  Th <- eval(parse(text=Th))
      Th <- as.double(Th)
      Th <- Th * 1e-9
    }
    #mus
    else if(  str_detect(Th, "(\\sμs)")  ) {
      Th <- str_remove_all(Th, "(\\sμs)")
      #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
      #  Th <- eval(parse(text=Th))
      Th <- as.double(Th)
      Th <- Th * 1E-6
    }
    #ms
    else if(  str_detect(Th, "(\\sms)")  ) {
      Th <- str_remove_all(Th, "(\\sms)")
      #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
      #  Th <- eval(parse(text=Th))
      Th <- as.double(Th)
      Th <- Th * 1e-3
    }
  #s
    else if(  str_detect(Th, "(\\ss)")  ) {
      Th <- str_remove_all(Th, "(\\ss)")
      #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
      #  Th <- eval(parse(text=Th))
      Th <- as.double(Th)
    }
  #min
  else if(  str_detect(Th, "(\\smin)")  ) {
    Th <- str_remove_all(Th, "(\\smin)")
    #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
    #  Th <- eval(parse(text=Th))
    Th <- as.double(Th) * 60
  }
    #hours
    else if(  str_detect(Th, "(\\sh)")  ) {
      Th <- str_remove_all(Th, "(\\sh)")
      #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
      #  Th <- eval(parse(text=Th))
      Th <- as.double(Th)
      Th <- Th * 60 * 60
    }
    #days
    else if(  str_detect(Th, "(\\sd)")  ) {
      Th <- str_remove_all(Th, "(\\sd)")
      #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
      #  Th <- eval(parse(text=Th))
      Th <- as.double(Th)
      Th <- Th * 60 * 60 * 24
    }
    #years
    else if(  str_detect(Th, "(\\sy)")  ) {
      Th <- str_remove_all(Th, "(\\sy)")
      #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
      #  Th <- eval(parse(text=Th))
      Th <- as.double(Th)
      Th <- Th * 60 * 60 * 24 * 365.25
    }
    #Th2 <- Th
    #rm(Th)
    return(Th)
}


##### MAKE THi and THf numeric #####
## Need to sort out the ~symbol at the beginning of an input
as_Type_Th_iso <- function(Thi) {
  if(is.na(Thi)) {return(NA)}
  ## account for the stable isomer
  if(  str_detect(Thi, "STABLE")  ){return("STABLE")}
  ## discount funny symbol cases
  # if(  str_detect(Thi, "@")  ){return(NA)}
  # #else {return(0)}
  # ## remove error and brackets
  # if( str_detect(Thi, "[<>~\\*\\+\\-]*") ) {
  #   Thi <- str_remove_all(Thi, "[<>~\\*\\+\\-]*")
  # }  
  # 
  # if(  str_detect(Thi, "(\\(\\d+\\))|[<>~\\*\\+]")  ){
  #   #if(str_detect(Thi, "<>~*")) {
  #   Th <- str_remove_all(Thi, "(\\(\\d+\\))|[<>~\\*\\+]") 
  #   if( str_detect(Th, "[<>~\\*\\+]")  ){
  #     #if(str_detect(Thi, "]<>~*]")) {
  #     Th <- str_remove_all(Th, "[<>~\\*\\+]") 
  #     Th <- trimws(Th)
  #   }  
  #   #}
  # }
  # 
  # else if( str_detect(Thi, "[<>~\\*\\+]")  ){
  #   #if(str_detect(Thi, "<>~*")) {
  #   Th <- str_remove_all(Thi, "[<>~\\*\\+]") 
  #   Th <- trimws(Th)
  # }  
  # else (Th <- Thi)
  Th <- Thi
  
  ## extract ns
  if(  str_detect(Th, "(\\sns)")  ) {
    Th <- str_extract_all(Th, "(\\sns)")

  }
  #mus
  else if(  str_detect(Th, "(\\sμs)")  ) {
    Th <- str_extract_all(Th, "(\\sμs)")
  }
  #ms
  else if(  str_detect(Th, "(\\sms)")  ) {
    Th <- str_extract_all(Th, "(\\sms)")
    
  }
  #s
  else if(  str_detect(Th, "(\\ss)")  ) {
    Th <- str_extract_all(Th, "(\\ss)")
  }
  #min
  else if(  str_detect(Th, "(\\smin)")  ) {
    Th <- str_extract_all(Th, "(\\smin)")
  }
  #hours
  else if(  str_detect(Th, "(\\sh)")  ) {
    Th <- str_extract_all(Th, "(\\sh)")
  }
  #days
  else if(  str_detect(Th, "(\\sd)")  ) {
    Th <- str_extract_all(Th, "(\\sd)")
  }
  #years
  else if(  str_detect(Th, "(\\sy)")  ) {
    Th <- str_extract_all(Th, "(\\sy)")
  }
  #Th2 <- Th
  #rm(Th)
  return(unlist(Th))
}




## Thi means any type of Th in this situation
as_Numeric_Th <- function(Thi) {
    if(is.na(Thi)) {return(NA)}

    ## assign stable states
    if(  str_detect(Thi, "STABLE")  ){return(-1)}

  ## convert "  ENERGY" to correct lifetime
  else if(  str_detect(Thi, "(\\sMEV)")  ) {
    Thi <- str_remove_all(Thi, "(\\sMEV)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi))
    Thi <- as.double(Thi)
    Thi <- hbar/(Thi*1E6*e)
  }
  ## convert "  ENERGY" to correct lifetime
  else if(  str_detect(Thi, "(\\sKEV)")  ) {
    Thi <- str_remove_all(Thi, "(\\sKEV)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi))
    Thi <- as.double(Thi)
    Thi <- hbar/(Thi*1E3*e)
  }
  ## convert "  ENERGY" to correct lifetime
  else if(  str_detect(Thi, "(\\sGEV)")  ) {
    Thi <- str_remove_all(Thi, "(\\sGEV)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi))
    Thi <- as.double(Thi)
    Thi <- hbar/(Thi*1E9*e)
  }
  ## convert "  ENERGY" to correct lifetime
  else if(  str_detect(Thi, "(\\sEV)")  ) {
    Thi <- str_remove_all(Thi, "(\\sEV)")
    #  Th <- paste(unlist(str_extract_all(Th, "[\\D+\\.]")), collapse = "")
    #  Th <- eval(parse(text=Th))
    Thi <- as.double(Thi)
    Thi <- hbar/(Thi*e)
   }
  
  
  #US
  else if(  str_detect(Thi, "(\\sUS)")  ) {
    Thi <- str_remove_all(Thi, "(\\sUS)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi)) 
    Thi <- as.double(Thi)
    Thi <- Thi * 1E-6
  }
  #FS
  else if(  str_detect(Thi, "(\\sFS)")  ) {
    Thi <- str_remove_all(Thi, "(\\sFS)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi)) 
    Thi <- as.double(Thi)
    Thi <- Thi * 1E-15
  }
  #PS
  else if(  str_detect(Thi, "(\\sPS)")  ) {
    Thi <- str_remove_all(Thi, "(\\sPS)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi)) 
    Thi <- as.double(Thi)
    Thi <- Thi * 1E-12
  }
  #NS
  else if(  str_detect(Thi, "(\\sNS)")  ) {
    Thi <- str_remove_all(Thi, "(\\sNS)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi)) 
    Thi <- as.double(Thi)
    Thi <- Thi * 1E-9
  }
      
    
    #ms
    else if(  str_detect(Thi, "(\\sMS)")  ) {
      Thi <- str_remove_all(Thi, "(\\sMS)")
      #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
      #  Thi <- eval(parse(text=Thi)) 
      Thi <- as.double(Thi)
      Thi <- Thi * 1e-3
    }
    #s
    else if(  str_detect(Thi, "(\\sS)")  ) {
      Thi <- str_remove_all(Thi, "(\\sS)")
      #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
      #  Thi <- eval(parse(text=Thi)) 
      Thi <- as.double(Thi)
    }
  #minutes
  else if(  str_detect(Thi, "(\\sM)")  ) {
    Thi <- str_remove_all(Thi, "(\\sM)")
    #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
    #  Thi <- eval(parse(text=Thi)) 
    Thi <- as.double(Thi) * 60
  }
    #hours
    else if(  str_detect(Thi, "(\\sH)")  ) {
      Thi <- str_remove_all(Thi, "(\\sH)")
      #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
      #  Thi <- eval(parse(text=Thi)) 
      Thi <- as.double(Thi)
      Thi <- Thi * 60 * 60
    }
    #days
    else if(  str_detect(Thi, "(\\sD)")  ) {
      Thi <- str_remove_all(Thi, "(\\sD)")
      #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
      #  Thi <- eval(parse(text=Thi)) 
      Thi <- as.double(Thi)
      Thi <- Thi * 60 * 60 * 24
    }
    #years
    else if(  str_detect(Thi, "(\\sY)")  ) {
      Thi <- str_remove_all(Thi, "(\\sY)")
      #  Thi <- paste(unlist(str_extract_all(Thi, "[\\D+\\.]")), collapse = "")
      #  Thi <- eval(parse(text=Thi)) 
      Thi <- as.double(Thi)
      Thi <- Thi * 60 * 60 * 24 * 365.25
    }
    return(Thi)
  }
##### MAKE THi and THf numeric #####  



  # ## nuclear photo-excitation search of a certain energy
# ## energy range available from photons energy "Exc"
# #Range <- 8 # keV
# #Exc <- 511
# #Qmin <- Exc - Range/2
# #Qmax <- Exc + Range/2
# 
# winExctf <- function(level, Einit, Qmin, Qmax) {
#   if(is.na(level$E) || is.na(Einit)) return(NA)
#   ## Evalutate nuclear excitation energy based on input
#   Q <- level$E - Einit
#   if(  (Qmin <= Q) && (Q <= Qmax)  ) { return(TRUE) }
#   else return(FALSE)
# }
# 

# 
# winExctfIso <- function(level, Einit, ind, Qmin, Qmax) {
#   if(is.na(level$E) || is.na(Einit)) return(NA)
#   ## Evalutate nuclear excitation energy based on input
#   Q <- level$E - Einit
#   if(  (Qmin <= Q) && (Q <= Qmax) ) { return(isomer$Key[ind]) }
#   else return(-1)
# }
# 
# 
# 



# ## TEST it returns level energy
# winBeamtfTEST <- function(level) {
#   return(level$E)
# }
# winBeamtfTEST(ND[4225,])

# ### TESTING ###
# Ax <- A0Ni_Kcapture$AX[1]
# Ax
# #vector of isomer indexes
# index <- which(isomer$AX==Ax)
# index
# winBeamtftup(A0Ni_Kcapture[1,], isomer$E[index[1]], index[1])
