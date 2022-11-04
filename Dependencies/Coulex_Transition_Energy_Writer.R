library(tidyverse)
library(stringr)

folder <-"~/w2k/Desktop/84Rb_Experiment/RbSiCoulex_CS_Estimate/"

file = paste0(folder,"84Rb_Levels.csv")
Levels <- read.csv(file, header=TRUE,fill=TRUE)

file <- paste0(folder,"84Rb_Transitions_gosia.csv")
Trans <- read.csv(file, header=TRUE, fill=TRUE)

Energies <- c(rep(0.0, length(Trans[,1])))
for(i in seq_along(Trans[,1])){
	Energies[i] <- Levels$Energy[Trans$Index2[i]] - Levels$Energy[Trans$Index1[i]]
}

sort(Energies)
