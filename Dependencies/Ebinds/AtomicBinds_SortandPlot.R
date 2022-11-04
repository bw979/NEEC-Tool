library(tidyverse)
library(stringr)
library(ggplot2)
library(plotly)

### Maybe need this 
## prototype ... discounts tentative assignment of Jpi
#ND_Eratio <- ND %>% filter(M > m1, M < m2, EvenEven == TRUE, Tentative_Jpi == FALSE) %>%
#  mutate(compute = integer(1), identify = logical(1), level1 = logical(1), level2 = logical(1), output = double(1)) 

#titles
FILE <- "Dependencies/Ebinds/NIST_Ebinds.csv"
Ebinds <- read_csv(FILE)

## Make a Jat column



### PLOTTING ###
#### EDIT HERE TO MAKE IT LOOK NICE ####
# p <- plot_ly(Ebinds, x =~Z, y =~N, z=~output, type="scatter3d", mode="markers+lines", line=list(width=2, dash="dash"), opacity=1.0, size=4, name=~El)
# htmlwidgets::saveWidget(p, title)



### Data Wrangling  ###
# ##ONCE: Convert data types
# #Ebinds$IE[1] <- trimws(Ebinds$IE[1], which = "left")
# Ebinds$IE <- as.double(Ebinds$IE)
# Ebinds$IE_uncertainty <- as.double(Ebinds$IE_uncertainty)
# Ebinds$Z <- as.double(Ebinds$Z)
# Ebinds$CS <- as.double(Ebinds$CS) 
#   
# ##ONCE: Rename headings and convert to keV
# Ebinds$IE <- Ebinds$IE *1E-3
# Ebinds$IE_uncertainty <- Ebinds$IE_uncertainty * 1E-3
# 
# write.table(Ebinds, file=FILE, append=T, row.names=F, col.names=T,  sep=",")


