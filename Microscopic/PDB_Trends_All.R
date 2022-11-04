library(tidyverse)
library(plotly)
#for the dfactorial function
library(phangorn)
library(stringr)
library(plotly)
library(XML)
library(processx)
library(pbapply)

#### conv for addition of Z ####
Conv <- read.csv("Dependencies/Z_Conversion.csv")
conv <- data.frame(Z=Conv[1], Symb=Conv[2], Name=Conv[3])
conv$Symb <- toupper(conv$Symb)

source("Wins.R")

# ###Calculated Data
wap_N0 <- read_csv("App_Data/All_RatesCalcd_Viva.csv", col_types = cols( Ji = col_character(), ICC = col_double(), ICC_tot = col_double(), S = col_double(), S_tot = col_double(), N_neec_pw_Stot = col_double(), N_neec_pw = col_double()   ))
wap <- wap_N0

## Consecutive Binding Energies
Ebinds <- read.csv("Dependencies/Ebinds/NIST_Ebinds.csv", header=TRUE, fill=TRUE)
## Neutral and consecutive binding energies with multipliers
ICC_Multipliers <- read_csv("Dependencies/ICC_Multipliers.csv")

### START HERE PRESENTLY
case_file <-"Microscopic/Theory_PDB_GS.csv"

#function declaration, no manipulation 
source("Microscopic/PDB_ResonanceStrength_Calculator.R")
source("Microscopic/WeisskopfCalculator.R")
source("Microscopic/ICC_Eval.R")

T_PDB <- read_csv(case_file)

### OCC is the number of electrons on the ion after capture
T_PDB$Occ <- T_PDB$Occ + 1
## Mutate to have a IC_SUMp column for gamma branching fraction
T_PDB <- mutate(T_PDB, ICC=double(1), ICC_tot=double(1), ICC_sum=double(1), Ei=0, Ef=double(1), Gamma_Scale_Factor=double(1) )
T_PDB <- select(T_PDB, AX, Z, El, A,Ei, Ef, Q, Vi, -Vi_q0, Ee, Ji, Jd, Jat, NEECX, Occ, -Div, Type, tau_SI, Gamma_eV, 
                oneoveroneplusalpha, -AroverGamma, B_wu, B_wu_2, Yneec_th_SI, Yneec_calc_SI, Ar_SI, Ar_AU, ICC, ICC_tot, ICC_sum, -IC_P, -IC_tot, -pi2_p2, -pi2_p2_th, 
                Sth, Ref, Spdb_P, -Spdb_tot, Ratio_P, -Rneec )




#### AnN0 and AnN0-X
## Calculate ICC's using up to date ICC calculator
for(i in 1:length(T_PDB$AX)){
   #T_PDB$ICC[i] <- filter(wap, AX==T_PDB$AX[i], Occ==T_PDB$Occ[i], Q==T_PDB$Q[i])$ICC
   T_PDB$ICC[i] <- Compute_ICC_singles(T_PDB$Type[i], T_PDB$Z[i], T_PDB$Q[i], (T_PDB$Occ[i]), "p", FALSE, 1, "M1+E2")
   T_PDB$ICC_tot[i] <- Compute_ICC_singles(T_PDB$Type[i], T_PDB$Z[i], T_PDB$Q[i], (T_PDB$Occ[i]), "tot", FALSE, 1, "M1+E2")
   T_PDB$ICC_sum[i] <- Compute_ICC_singles(T_PDB$Type[i], T_PDB$Z[i], T_PDB$Q[i], (T_PDB$Occ[i]), "SUMp", FALSE, 1, "M1+E2")
   T_PDB$Ar_SI[i] <- compute_Ar(T_PDB$Type[i], T_PDB$A[i], T_PDB$Q[i], T_PDB$B_wu[i], FALSE, 1)
}


## Calculate Resonance strength with gamma branching from SUMp already included
for(i in 1:length(T_PDB$AX)){
  T_PDB$Spdb_P[i] <- (1/(1+T_PDB$ICC_sum[i]))*compute_S(T_PDB$Jd[i],T_PDB$Jat[i],T_PDB$Ji[i],T_PDB$ICC[i], T_PDB$Ar_SI[i], T_PDB$Q[i], T_PDB$Vi[i])
  T_PDB$Yneec_calc_SI[i] <- compute_S_check_Y(T_PDB$Jd[i],T_PDB$Jat[i],T_PDB$Ji[i],T_PDB$ICC[i], T_PDB$Ar_SI[i], T_PDB$Q[i], T_PDB$Vi[i])
  T_PDB$Ratio_P[i] <- T_PDB$Spdb_P[i]/T_PDB$Sth[i]
}


# filename <- "Microscopic/Theory_PDB_2022_NEEC+X.csv"
# write.table(T_PDB, file=filename, append=T, row.names=F, col.names=T,  sep=",")  

#### AnNi