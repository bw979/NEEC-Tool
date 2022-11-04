library(tidyverse)
library(stringr)

## DATA SORTING
#NucLevels <- read.csv("~/w2k/Desktop/COULEX/RbSi2/84Rb_Levels0.csv", header=TRUE, fill=TRUE)
#NucLevels2 <- NucLevels %>% select(-Energy0)
#write.table(NucLevels2, file='~/Documents/current/COULEX/RbSi2/84Rb_Levels.csv', append=T, row.names=F, col.names=T,  sep=",")

#Trans <- read.csv("~/Documents/current/COULEX/RbSi2/84Rb_NNDC_Transitions2.csv", header=TRUE, fill=TRUE)



#Trans <- tibble(Index1=c(NA), Index2=c(NA), Energy=c(NA), deltaJ=c(NA), deltaPI=c(NA), Multipolarity=c(NA), ME=c(NA))

# Trans <- Trans %>% mutate(Energy = c(NA))
# Trans <- Trans %>% mutate(deltaJ = c(NA), deltaPI = c(NA), Multipolarity = c(NA), ME = c(NA))

## Set A
A <- 84

# Compute the Weisskopf estimates
BE1 <- 6.446E-4*A^(2/3)
BE2 <- 5.940E-6*A^(4/3)
BE3 <- 5.940E-8*A^2
BE4 <- 6.285E-10*A^(8/3)
BE5 <- 6.928E-2*A^(10/3)

### Where has this BMlam function come from?
BMlam <- function(lam) {(10*7/22)*(1.2)^(2*lam - 2)*((3/(lam + 3))^2)*A^((2*lam - 3)/3)}
#BM1 <- BMlam(1)
#BM1 <- 1.790
#BM2 <- BMlam(2)
#BM3 <- BMlam(3)
#BM4 <- BMlam(4)
BM1 <- 1.790
BM2 <- 1.650*A^(2/3)
BM3 <- 1.650*A^(4/3)
BM4 <- 1.746*A^2
BM5 <- 1.924*A^(8/3)

WE <- tibble(Type=character(9), one_Weisskopf_unit=double(9))
WE$Type[1] <- "BE1"
WE$Type[2] <- "BE2"
WE$Type[3] <- "BE3"
WE$Type[4] <- "BE4"
WE$Type[5] <- "BM1"
WE$Type[6] <- "BM2"
WE$Type[7] <- "BM3"
WE$Type[8] <- "BM4"
WE$Type[9] <- "BM5"

WE$one_Weisskopf_unit[1] <- BE1
WE$one_Weisskopf_unit[2] <- BE2
WE$one_Weisskopf_unit[3] <- BE3
WE$one_Weisskopf_unit[4] <- BE4
WE$one_Weisskopf_unit[5] <- BM1
WE$one_Weisskopf_unit[6] <- BM2
WE$one_Weisskopf_unit[7] <- BM3
WE$one_Weisskopf_unit[8] <- BM4
WE$one_Weisskopf_unit[9] <- BM5

write.table(WE, file='/mnt/H/Desktop/84Rb_Experiment/RbSiCoulex_CS_Estimate/Weisskopf_Estimates.csv', append=T, row.names=F, col.names=T,  sep=",")




#for(i in 1:length(NucLevels$Level)) {
  
  
#}


## FIll in the transition energies, deltaJ, deltaPI and subsequent Multipolarity
# for(i in 1:length(Trans$Index1)) {
#   if((Trans$Index2[i]) != 0){
#     Trans$Energy[i] = NucLevels$Energy[Trans$Index2[i]] - NucLevels$Energy[Trans$Index1[i]]
#     #change Parity?
#     if(NucLevels$Parity[Trans$Index2[i]] == NucLevels$Parity[Trans$Index1[i]]) {Trans$deltaPI[i] = FALSE} else {Trans$deltaPI[i] = TRUE}
#     #deltaJ
#     Trans$deltaJ[i] = abs(NucLevels$Spin[Trans$Index2[i]] - NucLevels$Spin[Trans$Index1[i]])
#     #Fills in multipolarity
#     switch(Trans$deltaJ[i],
#            if(Trans$deltaPI[i] == TRUE) Trans$Multipolarity[i] = 'E1' else Trans$Multipolarity[i] = 'M1',
#            if(Trans$deltaPI[i] == TRUE) Trans$Multipolarity[i] = 'M2' else Trans$Multipolarity[i] = 'E2',
#            if(Trans$deltaPI[i] == TRUE) Trans$Multipolarity[i] = 'E3' else Trans$Multipolarity[i] = 'M3',
#            if(Trans$deltaPI[i] == TRUE) Trans$Multipolarity[i] = 'M4' else Trans$Multipolarity[i] = 'E4')
#      # if(Trans$Multipolarity[i] == 'NA') Trans$ME = NA
#      
#      # else Trans$ME = -1
#     }
# }


for(i in 1:length(Trans$Index1)) {
  if((Trans$Index2[i]) != 0) {
    if(is.na(Trans$Multipolarity[i])) Trans$ME[i] = 0
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'E1') Trans$ME[i] = BE1
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'E2') Trans$ME[i] = BE2
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'E3') Trans$ME[i] = BE3
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'E4') Trans$ME[i] = BE4
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'M1') Trans$ME[i] = BM1
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'M2') Trans$ME[i] = 0.0
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'M3') Trans$ME[i] = 0.0
    if(!is.na(Trans$Multipolarity[i]) & Trans$Multipolarity[i] == 'M4') Trans$ME[i] = 0.0
  }
}
 

#Trans2 <- Trans %>% filter(deltaJ != 0, ME != -1)
#Trans2 <- Trans %>% filter(deltaJ != 0 )

Trans2 <- Trans %>% select(Index1, Index2, ME)

#if(Trans2$ME[i] = NA) {}


## Write the Trans to file
write.table(Trans2, file='84Rb_MatrixElements.csv', append=T, row.names=F, col.names=F,  sep=",")






# ## Add the deltaJ column
# for(i in 1:length(Trans$Index1)){
#   Trans$deltaJ[i] = NucLevels$Spin[Trans$deltaJ[i]]
# }


# for(i in 1:length(Trans$Index1)) {
#   if((Trans$Index2[i]) != 0){
#     switch(Transitions$deltaJ[i],
#       if(Trans$deltaPI[i] == 2) Trans$Multipolarity[i] = 'E1' else Trans$Multipolarity[i] = 'M1',
#       if(Trans$deltaPI[i] == 2) Trans$Multipolarity[i] = 'M2' else Trans$Multipolarity[i] = 'E2',
#       if(Trans$deltaPI[i] == 2) Trans$Multipolarity[i] = 'E3' else Trans$Multipolarity[i] = 'M3',
#       if(Trans$deltaPI[i] == 2) Trans$Multipolarity[i] = 'M4' else Trans$Multipolarity[i] = 'E4'
#     )
#   }
# }




