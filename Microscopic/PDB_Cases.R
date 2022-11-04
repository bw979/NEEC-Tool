############################# CASES #################################
#Need to modify these, read in from the EBind list, then do a final Palffy comparison 

##84Rb isomer case
Eres_eV <- 3050
BW <-0.08
Ji <-6
Jd <-5
Jat <-0.5
IC <- 0.32
Eres_MeV <- Eres_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- 1.790 * BW
Ar <- 1.779E13*(Eres_MeV^3)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,Eres_eV)
### OUTPUT : 0.01066692 ... Agrees with same calc done by Palffy


#### SOME PALFFY HYDROGENIC CASES 
#165Ho case
En_eV <- 94700
BW <- 0.275
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- 1.790 * BW
Ar <- 1.779E13*(En_MeV^3)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
IC <- 3.0
compute_S(7/2,0.5,9/2,IC,Ar,En_eV)
### OUTPUT : 10.09123 beV ... Palffy Value:0.884 beV


# ## NO B or IC VALUE #### #173Yb case
# Eres_eV <- 78640
# BW <- 0.275
# SP_WeisskopfM1_SI <- 1.79
# Eres_MeV <- Eres_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- 1.790 * BW
# Ar <- 1.779E13*(Eres_MeV^3)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# IC <- 3.0
# compute_S(7/2,0.5,5/2,IC,Ar,Eres_eV)
# ### OUTPUT : 

# ## NO B or IC VALUE #### #185Re case
# Eres_eV <- 
# BW <- 0.275
# SP_WeisskopfM1_SI <- 1.79
# Eres_MeV <- Eres_eV*(10^-6)
# #Ar_TM1 <- 1.779E13*
# B_SI <- 1.790 * BW
# Ar <- 1.779E13*(Eres_MeV^3)*B_SI # PER SECOND ... works spot on
# #compute_Ar(1,B_SI,3200) #  Should match
# #Ar_SI <- Ar / time_conv
# IC <- 3.0
# compute_S(7/2,0.5,5/2,IC,Ar,Eres_eV)
# ### OUTPUT : 

#187Re case
En_eV <- 134000.243
BW <-0.260
Ji <-5/2
Jd <-7/2
Jat <-0.5
IC <- 2.2
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- 1.790 * BW
Ar <- 1.779E13*(En_MeV^3)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
### OUTPUT : 16 beV ... Palffy Value:1.16 beV



#####  55Mn M1
En_eV <- 125.95E3
BW <-0.0417
Ji <-5/2
Jd <-7/2
Jat <-0.5
IC <- 0.01691
Eres_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- 1.790 * BW
Ar <- 1.779E13*(En_MeV^3)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
### L shell IC OUTPUT : 0.01903 beV ... Palffy Value: 9.22E-4 beV
### Total IC OUTPUT :   0.01903 beV




##### E2 #####

#####  129Sb 17min isomer E2
En_eV <- 9.76E3
BW <-1.96
Ji <-19/2
Jd <-15/2
Jat <-1.5
IC <- 5.59E3
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 129^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)
### L shell OUTPUT : 0.01903 beV ... Palffy Value: 1.6E-3 beV
### M shell OUTPUT : 0.00586 beV ... Palffy Value: 3E-5 beV


## 164 Dy
En_eV <- 73.392E3
BW <- 211
Ji <-0
Jd <-2
Jat <-0.5
IC <- 8.89
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 164^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)


## 170Er
En_eV <- 78.591E3
BW <- 208
Ji <-0
Jd <-2
Jat <-0.5
IC <- 7.47
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 164^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)

## 174Yb
En_eV <- 76.471E3
BW <- 201
Ji <-0
Jd <-2
Jat <-0.5
IC <- 9.43
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 164^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)


## 174Yb
En_eV <- 76.471E3
BW <- 201
Ji <-0
Jd <-2
Jat <-0.5
IC <- 9.43
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 164^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)

## 154Gd
En_eV <- 123.071E3
BW <- 157
Ji <-0
Jd <-2
Jat <-0.5
IC <- 2.297
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 164^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)

## 156Gd
En_eV <- 88.966E3
BW <- 189
Ji <-0
Jd <-2
Jat <-1.5
IC <- 0.24
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 164^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)

## 62Dy
En_eV <- 80.066E3
BW <- 204
Ji <-0
Jd <-2
Jat <-0.5
IC <- 6.14
En_MeV <- En_eV*(10^-6)
#Ar_TM1 <- 1.779E13*
B_SI <- (5.94E-2)* 164^(4/3) * BW
Ar <-(1.223E9)*(En_MeV^5)*B_SI # PER SECOND ... works spot on
#compute_Ar(1,B_SI,3200) #  Should match
#Ar_SI <- Ar / time_conv
compute_S(Jd,Jat,Ji,IC,Ar,En_eV)










# #Data to plot
# df=ebinds
# 
# #### EDIT HERE TO MAKE IT LOOK NICE ####
# p <- plot_ly(df, x =~Z, y =~, z=~output, type="scatter3d", mode="markers+lines", line=list(width=2, dash="dash"), opacity=1.0, size=4, name=~El)
# htmlwidgets::saveWidget(p, "EvenEven_Plots.html")




## Try a few other known cases
##












# c <- 137
# Etrans <- 3000 #eV
# Etrans_J <- Etrans_SI_keV *e
# Etrans_AU <- Etrans_SI_J/Eh
# ##other way
# E_AU <- 110
# E_J <- E_AU * Eh
# E_eV <- E_J / e
# ## ... works
# p <- Etrans_AU / c
# Spdb_AU <- compute_S(5,0.5,6,3400,Ar)





# The p bit
# m <- 1
# E <- p/c in atomic units 
# # stay in atomic units then one big conversion at the end



### ... ideally
#A0N0_S_pdb <- compute_S(A0N0_S_pdb)

##compare these to palffys theoretical estimates in her thesis using a well reasoned fraction


#filter these observations where Ed - Ei = Egam => keep IC    .... A0N0_S_GAM$Etrans <- round(A0N0_S_GAM$Ed - A0N0_S_GAM$Ei, 2)

# ### Current semi shit way of doing it
# for(i in 1:length(A0N0_S_GAM$Etrans[i])){
#   #if the transition matches assign A0N0_S$IC[i] <- CORRECT ND_Gamma$CC[i]
#   
#   Etrans <- round((A0N0_S$Ed[i] - A0N0_S$Ei[i]), 2)
#   if( Etrans == ND_Gamma$E[i] ) {
#     A0N0_S$IC[i] <- A0N0_S$CC[i]
#   }
#   else {
#     A0N0_S$IC[i] <- NA
#   }
#   
# }
# 
# A0N0_S <- filter(!is.na(A0N0_S$IC))




# #### Returns a value TRUE of FALSE if the IC coeff is known for the NEEC transition
# gammaTrans <- function(Etrans, Egam) {
#   if( Etrans == Egam ){
#     return(TRUE)
#   }
#   else {
#     return(FALSE)
#   }
# }
# 
# A0N0_S_GAM$IC_NEEC <- map2(A0N0_S_GAM$Etrans, A0N0_S_GAM$Egam, gammaTrans)
# 






#### Taken from AnNi, finds and produces a vector of indexes
# for(i in 1:length(ND_AnNi$E)){
#   #store current isomer to be looked at
#   #progress(i)
#   Ax <- ND_AnNi$AX[i]
#   
#   #vector of isomer indexes
#   index <- which(isomer$AX==Ax)
# 
# }


###Assign IC column 
##Need to do another join here
#A0N0_S$IC <- filter(ND_Gamma, AX = ax)
# # try for loop
# for( i in 1:length(A0N0_S$IC) ) {
#   #if the transition matches assign A0N0_S$IC[i] <- CORRECT ND_Gamma$CC[i]
#   
#   #could use arr <-  which(ND_Gamma$AX == "12B")
#   #outputs : arr <- c(2 3 4)
#   
#   
#   #...shit way of doing it, nested for loops
#   #round(x, digits = 0)
#   for( j in 1:length(ND_Gamma)){
#     Etrans <- round((A0N0_S$Ed[i] - A0N0_S$Ei[i]), 2)
#     if( Etrans == ND_Gamma$E[j] ) {
#       A0N0_S$IC[i] <- ND_Gamma$CC[j]  
#     }
#   }
# }

#ax <- levels(as.factor(ND_Gamma$AX))

