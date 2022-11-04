n_index <- read_csv("Dependencies/n_to_ICshell_Conversion.csv")
ICC_Multipliers <- read_csv("Dependencies/ICC_Multipliers.csv")

############ Compute_ICC_sinlges is most up to date function #########
#### Feb 2022  - updated to read in subshell to ensure correct CCp is read out in K forbidden transitions

# i<-80
# Type <- Wap$Type[i]; Zin <- Wap$Z[i]; Q <- Wap$Q[i]; n <- (Wap$Occ[i]); OUTPUT <- "tot"
# Compute_ICC_singles(Wap$Type[i], Wap$Z[i], Wap$Q[i], 2, "tot")

#Wap <- filter(Wap, AX=="93MO", Q<4.86)
Compute_ICC_singles <- function(Type, ss,  Zin, Q, n, OUTPUT, MIX, MR, MIX_type){
   # i <- 1
   #Type <- Wap$Type[i]; Zin <- Wap$Z[i];  Q <- Wap$Q[i] ; n <- Wap$Occ[i]; OUTPUT <- "p"; MIX<-FALSE; MR <- Wap$MR[i]; MIX_type <- Wap$MIX_type[i]
  #Type <- Gamma$Type[i], ss<- "", Gamma$Z[i], Gamma$Egam[i], 1, "tot", FALSE, 1, ""
  #e_file <- sprintf("Dependencies/Ebinds/Configuration_Files/Z_%003.f.csv", Z)
  # i<-1
  # Type<-"E2"; Zin<-42; Q<-4.85; n<-Abinit$Occ_; OUTPUT<-"p"; MIX<-FALSE, MR<-0, MIX_type<-""; ss<-"M1"
  # 
  #ch <- x
  
  atom <- filter(ICC_Multipliers, Z==Zin)
  #if(Z == c(10, 15, 20, 25,30, 35, 40,50,60,70,80,90,100)) print(Z)
  if(is.na(Type) || Type == "Large_L" || Type == "E0" || Zin < 11 || Zin > 103)
  {
    return(NA)
  }
  else {
    if(MIX==FALSE){
      bricc_string <- sprintf("./briccs -Z %d -g %.5f -L %s -a -w BrIccFO >> output.xml", Zin, Q, Type)
      unlink("output.xml")
      system(bricc_string)
    } else {
      bricc_string <- sprintf("./briccs -Z %d -g %.5f -L %s -a -w BrIccFO >> output.xml", Zin, Q, Type)
      if(MIX_type == "M1+E2"){
        bricc_string <- sprintf("./briccs -Z %d -g %.5f -L %s -d %.5f -a -w BrIccFO >> output.xml", Zin, Q, MIX_type, MR)
      } else {bricc_string <- sprintf("./briccs -Z %d -g %.5f -L %s -a -w BrIccFO >> output.xml", Zin, Q, Type)}
      
      unlink("output.xml")
      system(bricc_string)
    }
    
    ## Parse the xml file to data.frame here
    ## Want to grab relevant IC coeff from n
    #n is 1 for now so grab IC_k
    #### ELEMENT 5 is IC_tot so IC_k is element 5 + 1 
    #### IC L2 would
    
    ###### OK Actually gonna have to read out all shells into a tibble then allocate
    data <- xmlParse("output.xml")
    x <- xmlToList(data)
    len <- length(x) - 6
    
    # if( unlist(x[5]$WARNING$text) == "  Transition energy is within 1 keV to one of the shell binding energies. " ){
    #   return(NA) 
    # }
    
    ## [6] is always the location of the first CC coeff
    ## TODO:
    # -n to shell conversion table as well
    # -parse the input string for briccs
    
    if(len > 0){
      
      Shells <- tibble(subshell = character(len), CCp = double(len), Eic = double(len) )
      ####MIXING FALSE###  
      if(MIX==FALSE){
        for(k in 1:len){
          #Shells$subshell[k] <- x[5+k]$MixedCC$.attrs[[1]]
          Shells$subshell[k]  <- x[5 + k]$PureCC$.attrs[[1]]
          Shells$Eic[k]  <- as.double(x[5 + k]$PureCC$.attrs[[2]])
          
          
          #else {
          ## Read the partial ICC into that shell in  the temporary Shells tibble
          ###### THIS FINALLY WORKS BECAUSE I unlisted the result of str_extract_all otherwise was producing a list containing 1 result
          ICCp <- paste(unlist(str_extract_all(x[5 + k]$PureCC$text, "[\\d\\.E-]")), collapse = "")
          ICCp <- eval(parse(text=ICCp))
          if( is.null(ICCp) ) {
            Shells$CCp[k] <- NA
            next
          }
          if(!is.null(ICCp)) {
            Shells$CCp[k] <- ICCp
          }
          
          #ICCtot <- paste(unlist(str_extract_all(x[5 + k]$PureCC$text, "[\\d\\.E-]")), collapse = "")
        }
        #close mix conditional
      }
      
      ####MIXING TRUE### 
      if(MIX==TRUE){
        for(k in 1:len){
          Shells$subshell[k] <- x[5+k]$MixedCC$.attrs[[1]]
          #Shells$subshell[k]  <- x[5 + k]$PureCC$.attrs[[1]]
          
          Shells$Eic[k]  <- as.double(x[5 + k]$MixedCC$.attrs[[2]])
          
          #else {
          ## Read the partial ICC into that shell in  the temporary Shells tibble
          ###### THIS FINALLY WORKS BECAUSE I unlisted the result of str_extract_all otherwise was producing a list containing 1 result
          ICCp <- paste(unlist(str_extract_all(x[5 + k]$MixedCC$text, "[\\d\\.E-]")), collapse = "")
          ICCp <- eval(parse(text=ICCp))
          if( is.null(ICCp) ) {
            Shells$CCp[k] <- NA
            next
          }
          if(!is.null(ICCp)) {
            Shells$CCp[k] <- ICCp
          }
          
          #ICCtot <- paste(unlist(str_extract_all(x[5 + k]$PureCC$text, "[\\d\\.E-]")), collapse = "")
        }
      } 
      
      atom <- filter(atom, subshell %in% Shells$subshell )
      
      #DCC <- x[6 + k]$PureCC$.attrs[1]
      
      if(OUTPUT == "p"){
        ### Now assign the required CC based on n as input
        ##### REMEMBER WE ARE TAKING THE AVERAGE SHELL VALUE #####
        #Check that ICC has been considered for that shell by briccs
        if(  atom$subshell[n] %in% Shells$subshell  ){
          
          ################ ss addition
           Occ <- n
           Multiplier <- atom$multiplier[which(atom$subshell == ss & atom$Occ == n )]
           ICCp <- Shells$CCp[which(Shells$subshell ==  ss)] #* Multiplier  
          
          ######################
          
       ##### USUAL:   
          # Occ <- n
          # Multiplier <- atom$multiplier[which(atom$Occ == n )]
          # ICCp <- Shells$CCp[which(Shells$subshell ==  atom$subshell[Occ])]# * Multiplier
        ###########   
          #ICCp <- Shells$CCp[which(Shells$subshell ==  n_index$subshell[n])] * Multiplier
        }
        else {
          ICCp <- NA
        }
        
        #Shell <- x[6]$PureCC$.attrs[[1]]
        #IC_value <- as.double(str_extract(x[6]$PureCC$text, "\\d+"))
        ## Extract the numeric value from the text string
        if( is.null(ICCp) || is.na(ICCp) ) {
          return(NA)
          #next
        }
        #if( is.list())
        #if( is.na(ICCp) ) return(NA)
        else {
          return(ICCp)
        }
      }
      
      if(OUTPUT == "FRAC"){
        #use atom to filter Shells between Occ=Occmin and Occmax-1
        #ICCfrac_SUM <- 0
        
        ##TODO
        ### Need to add nlj subshell averaging fraction ... DONT NEED BY COMPARING TO NEECx PAPER ... BUT SHOULD ACTUALLY IMPLEMENT ANOTHER WEIGHT FUNCTION
        
        ### Need to add starting index if it is K forbidden ... ie the first index in ICC_multipliers that is from Shells$subshel
        min_n <- atom$Occ[which(atom$subshell==Shells$subshell[1])]
        
        for(k in min_n:n){
          #k<-80
          #if(  n_index$subshell[n] %in% Shells$subshell  ){
          if(  atom$subshell[k] %in% Shells$subshell  ){
            ## Sums only the subshells that are in the atom list filtered below occupation
            ICCfrac_SUM <- sum(filter(Shells, subshell %in%  filter(atom, Occ<=n)$subshell)$CCp)
            #ICCp <- Shells$CCp[which(Shells$subshell ==  n_index$subshell[n])] * Multiplier
          }
          else {
            #ICCp <- NA
            ICCfrac_SUM <- NA
          }
        }
        
        #Shell <- x[6]$PureCC$.attrs[[1]]
        #IC_value <- as.double(str_extract(x[6]$PureCC$text, "\\d+"))
        ## Extract the numeric value from the text string
        if( is.null(ICCfrac_SUM) || is.na(ICCfrac_SUM) ) {
          return(NA)
          #next
        } else {
          return(ICCfrac_SUM)
        }
        
      }
      
      if(OUTPUT == "SUMp"){
        ### Now assign the required CC based on n as input
        ##### REMEMBER WE ARE TAKING THE AVERAGE SHELL VALUE #####
        #Check that ICC has been considered for that shell by briccs
        ICCpSUM <- 0
        
        ##TODO
        ### Need to add nlj subshell averaging fraction ... DONT NEED BY COMPARING TO NEECx PAPER ... BUT SHOULD ACTUALLY IMPLEMENT ANOTHER WEIGHT FUNCTION
        
        ### Need to add starting index if it is K forbidden ... ie the first index in ICC_multipliers that is from Shells$subshel
        min_n <- atom$Occ[which(atom$subshell==Shells$subshell[1])]
        
        for(k in min_n:n){
          #k<-80
          #if(  n_index$subshell[n] %in% Shells$subshell  ){
          if(  atom$subshell[k] %in% Shells$subshell  ){
            #Occ <- k
            ### This is is the multiplier for that subshell times the averaging of the nl over that nlj
            Multiplier <- atom$multiplier[which(atom$Occ == k )]#*(1/atom$nmax[which(atom$Occ == k )])
            #nnmax <- atom$nmax[which(atom$Occ == n )]
            ICCp <- Shells$CCp[which(Shells$subshell ==  atom$subshell[k])] * Multiplier
            ICCpSUM <- ICCpSUM + ICCp
            #ICCp <- Shells$CCp[which(Shells$subshell ==  n_index$subshell[n])] * Multiplier
          }
          else {
            #ICCp <- NA
            ICCpSUM <- NA
          }
        }
        
        #Shell <- x[6]$PureCC$.attrs[[1]]
        #IC_value <- as.double(str_extract(x[6]$PureCC$text, "\\d+"))
        ## Extract the numeric value from the text string
        if( is.null(ICCpSUM) || is.na(ICCpSUM) ) {
          return(NA)
          #next
        } else {
          return(ICCpSUM)
        }
      }
      
      
    
      
      #}
      if(OUTPUT == "tot"){
        ###### BIT FOR ADDING ICCtot
        #if(len > 0){
        if(MIX=="FALSE"){ICCtot <- paste(unlist(str_extract_all(x[5]$PureCC$text, "[\\d\\.E-]")), collapse = "")}
        if(MIX=="TRUE"){ICCtot <- Shells$CCp[1]}
        ICCtot <- eval(parse(text=ICCtot))
        if( is.null(ICCtot) ) {
          return(NA)
          #next
        }
        if(!is.null(ICCtot)) ICC <- ICCtot
        
        #if(  is.null(ICC)  ) {
        #  return(NA)
        #  #next
        #}
        #if( is.list())
        if( is.na(ICC) ) return(NA)
        else {return(ICC)}
        #}
        #else {return(NA)} 
      }
    }
    else {
      return(NA)
    }
    
        
      
    # ### Return The neutral binding energy of that subshell
    # if(OUTPUT == "Eic"){
    #   if(  n_index$subshell[n] %in% Shells$subshell  ){
    #     Eic <- Shells$Eic[n]
    #   }
    #   else {
    #     Eic <- NA
    #   }
    # 
    #   #Shell <- x[6]$PureCC$.attrs[[1]]
    #   #IC_value <- as.double(str_extract(x[6]$PureCC$text, "\\d+"))
    #   ## Extract the numeric value from the text string
    #   if( is.null(Eic) || is.na(Eic) ) {
    #     return(NA)
    #     #next
    #   }
    #   #if( is.list())
    #   #if( is.na(ICCp) ) return(NA)
    #   else {
    #     return(Eic)
    #   }
    #  }
    
    # }
    #return()
    
  }
  
  ##End not mixing statement  
  #}  
  
  
  ##### END FUNCTION    
}

      
