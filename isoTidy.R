library(tidyverse)
library(stringr)

unlink("isomer.csv")
rm(isomer)

## TODO: # Get rid of funny symbol in isomer AX

#AtNI
AtNI <- read_csv("Dependencies/AtNI.csv")
AtNI$AX <- toupper(AtNI$AX)
#ISOMERS
isomer <- data.frame(AX=AtNI$AX, M=(as.numeric(as.character(AtNI$Z)) + as.numeric(as.character(AtNI$N))), Z=AtNI[1], E=AtNI[4], J=AtNI[5], Thalf=AtNI[6], Egam=AtNI[7], lam=AtNI[8])
isomer <- transform(isomer, E.keV.=gsub("\\s[0-9()]+","", E.keV.))
isomer <- filter(isomer, !is.na(M))
isomer <- as_tibble(isomer)
isomer$E.keV. <- as.double(isomer$E.keV.) 

isomer$Z <- as.integer(isomer$Z)
isomer$M <- as.integer(isomer$M)

isomer$AX <- as.character(isomer$AX)
isomer <- isomer %>% rename(E = E.keV.)

#T col
isomer <- isomer %>% mutate(TF = TRUE)

#Key Col
isomer <- isomer %>% mutate(Key = " ")
isomer$Key <- with(isomer, paste0(E, ";", J, ";", T1.2))

isomer <- isomer %>% filter(!is.na(E))
isomer <- isomer %>% rename(Thalf = T1.2)



#write CSV
if(!file.exists("Dependencies/isomer.csv")){
write.table(isomer, file="Dependencies/isomer.csv", append=T, row.names=F, col.names=T,  sep=",")
}

