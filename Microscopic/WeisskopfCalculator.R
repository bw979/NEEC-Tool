###### Weisskopf Ar calculator #####
### Compute the Ar column based on tabulated B values
#WC <- read.csv("WE.csv")
#Evaluate AR_AU using WE table
## THis will be a function that outputs the Ar when given a B value in WU along with the relevant A, Etrans and type



#Etrans in keV
compute_Ar <- function(Type, A, Etrans, B_wu, MIX, MR) {
  if(is.na(Type)) return(NA)
  if(Type == "E1"){
    Ar_SI <- (1.587E15) * ((Etrans * 1E-3)^3) * (6.446E-2)*((A)^(2/3)) * B_wu
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "E2"){
    Ar_SI <- (1.223E9) * ((Etrans * 1E-3)^5) * (5.940E-2)*((A)^(4/3)) * B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "E3"){
    Ar_SI <- (5.698E2) * ((Etrans * 1E-3)^7) *  (5.940E-2)*((A)^(2)) * B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "E4"){
    Ar_SI<- (1.694E-4) * ((Etrans * 1E-3)^9) * (6.285E-2)*((A)^(8/3)) * B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "E5"){
    Ar_SI <- (3.451E-11) * ((Etrans * 1E-3)^11) * (6.928E-2)*((A)^(10/3)) *B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "M1"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    Ar_SI <- (1.779E13) * ((Etrans * 1E-3)^3) * 1.79* B_wu * (A/A) * (B_wu/B_wu)
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "M2"){
    Ar_SI <- (1.371E7) * ((Etrans * 1E-3)^5) * 1.650*((A)^(2/3)) *B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "M3"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    Ar_SI <- (6.387E0) * ((Etrans * 1E-3)^7) *  1.650*((A)^(4/3)) *B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "M4"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    Ar_SI <- (1.899E-6) * ((Etrans * 1E-3)^9) * 1.746*((A)^(2)) *B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  if(Type == "M5"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    Ar_SI <- (3.868E-13) * ((Etrans * 1E-3)^11) * 1.924*((A)^(8/3)) *B_wu 
    Ar_AU <- Ar_SI * time_conv
  }
  
  if(MIX==FALSE){
    return(Ar_SI)
  }
  if(MIX==TRUE){
    return(Ar_SI + ((MR^2)*Ar_SI) )
  }
}

## computes the matrix elements in units e.b^lambda/2 from an input B(lambda) in weisskopf units
compute_WE_eb <-  function(Type, A, B_wu, Ii, Id) {
  if(is.na(Type)) return(NA)
  if(Type == "E1"){
    ## e^2 b^lambda
    B_down <- (6.446E-4)*((A)^(2/3)) * B_wu
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "E2"){
    B_down <- (5.940E-6)*((A)^(4/3)) * B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "E3"){
    B_down <-  (5.940E-8)*((A)^(2)) * B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "E4"){
    B_down <- (6.285E-10)*((A)^(8/3)) * B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "E5"){
    B_down <-  (6.928E-12)*((A)^(10/3)) *B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "M1"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    B_down <-  1.79* B_wu * (A/A) * (B_wu/B_wu)
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "M2"){
    B_down <- 1.650E-2*((A)^(2/3)) * B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "M3"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    B_down <-  1.650E-4*((A)^(4/3)) *B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "M4"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    B_down <-  1.746E-6*((A)^(2)) *B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  if(Type == "M5"){
    #Ar_AU <- WC$Ar_AU[ which(WC$Type=="M1") ]
    B_down <- 1.924E-8*((A)^(8/3)) *B_wu 
    B_up <- ((2*Id + 1)/(2*Ii + 1)) * B_down
    ME <- sqrt((2*Ii + 1)*B_up)
  }
  return(ME)
}
  

  
  

  

# compute_Ar_MIXED <- function(Type, A, Etrans, B_wu, B_wu2, delta) {
# #B1
#   Ar_B1 <- compute_Ar(Type, A, Etrans, B_wu)
#   
# #B2  
#   Ar_B2 <- compute_Ar(Type, A, Etrans, B_wu2)
#   
#   
#   Ar_TOT <-
# }









###### Weisskopf Ar calculator #####
### Compute the Ar column based on tabulated B values
#WC <- read.csv("WE.csv")
#Evaluate AR_AU using WE table
## THis will be a function that outputs the Ar when given a B value in WU along with the relevant A, Etrans and type


# for(i in 1: length(T_PDB$Type)){ 
# print(i)
# if(T_PDB$compute[i] == TRUE){
#   if(T_PDB$Type[i] == "E1"){
#     T_PDB$Ar_SI[i] <- (1.587E15) * ((T_PDB$Etrans[i] * 1E-3)^3) *      (6.446E-4)*((T_PDB$A[i])^(2/3)) * T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "E2"){
#     T_PDB$Ar_SI[i] <- (1.223E9) * ((T_PDB$Etrans[i] * 1E-3)^5) *     (5.940E-6)*((T_PDB$A[i])^(4/3)) * T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "E3"){
#     T_PDB$Ar_SI[i] <- (5.698E2) * ((T_PDB$Etrans[i] * 1E-3)^7) *         (5.940E-8)*((T_PDB$A[i])^(2)) * T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "E4"){
#     T_PDB$Ar_SI[i] <- (1.694E-4) * ((T_PDB$Etrans[i] * 1E-3)^9) *       (6.285E-10)*((T_PDB$A[i])^(8/3)) * T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "E5"){
#     T_PDB$Ar_SI[i] <- (3.451E-11) * ((T_PDB$Etrans[i] * 1E-3)^11) *       (6.928E-2)*((T_PDB$A[i])^(10/3)) * T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "M1"){
#     #T_PDB$Ar_AU[i] <- WC$Ar_AU[ which(WC$Type=="M1") ]
#     T_PDB$Ar_SI[i] <- (1.779E13) * ((T_PDB$Etrans[i] * 1E-3)^3) *          1.79* T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "M2"){
#     T_PDB$Ar_SI[i] <- (1.371E7) * ((T_PDB$Etrans[i] * 1E-3)^5) *        1.650*((T_PDB$A[i])^(2/3)) *T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "M3"){
#     #T_PDB$Ar_AU[i] <- WC$Ar_AU[ which(WC$Type=="M1") ]
#     T_PDB$Ar_SI[i] <- (6.387E0) * ((T_PDB$Etrans[i] * 1E-3)^7) *       1.650*((T_PDB$A[i])^(4/3)) *T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   
#   if(T_PDB$Type[i] == "M4"){
#     #T_PDB$Ar_AU[i] <- WC$Ar_AU[ which(WC$Type=="M1") ]
#     T_PDB$Ar_SI[i] <- (1.899E-6) * ((T_PDB$Etrans[i] * 1E-3)^9) *         1.746*((T_PDB$A[i])^(2)) *T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
#   if(T_PDB$Type[i] == "M5"){
#     #T_PDB$Ar_AU[i] <- WC$Ar_AU[ which(WC$Type=="M1") ]
#     T_PDB$Ar_SI[i] <- (3.868E-13) * ((T_PDB$Etrans[i] * 1E-3)^11) *        1.924*((T_PDB$A[i])^(8/3)) *T_PDB$B_wu[i] 
#     T_PDB$Ar_AU[i] <- T_PDB$Ar_SI[i] * time_conv
#   }
# }
# }
#WE_Ar <- read_csv("Dependencies/WE.csv")
##Ar - Radiative transition rate (first read B(lam) from GAMMA continuation card)
#All inputs in AU 
# compute_Ar <- function(type, B, Etrans){
#   #Etrans_MeV <- Etrans_eV*(10^-6)
#   Ar <- WE_Ar$Ar[which(WE_Ar$Type==type)]
#   Ar_AU <- Ar / time_conv
#   return(Ar_AU) 
# }

#RETURNS W-Eckart estimate function (SI units atm)
#INPUT: L, B(multiples of W.U.), Etrans (eV)
# compute_Ar <- function(L, B, Etrans){
#   #Etrans_MeV <- Etrans_eV*(10^-6)
#   c <- 137.036
#   Ar <- ((8*pi*(L+1))/(L*((dfactorial(2*L+1))^2)))*((Etrans^(2*L+1))/(c))*B
#   return(Ar)
# }







