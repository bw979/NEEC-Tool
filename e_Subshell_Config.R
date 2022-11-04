library(RColorBrewer)
source("Wins.R")
n_index <- read_csv("Dependencies/n_to_ICshell_Conversion.csv")
Atomic_Colours <- read_csv("Dependencies/Atomic_Colours.csv")
#### conv for identifying an element from Z ####
Conv <- read.csv("Dependencies/Z_Conversion.csv")
conv <- data.frame(Z=Conv[1], Symb=Conv[2], Name=Conv[3])
conv$Symb <- toupper(conv$Symb)
whi <- function(a) {return(which(conv$Symb==a))}


## Look at the consecutive ionisation energies
#Z_choose <- 37

for(Z_choose in 10:103){

## Unlink the file in case its a rewrite
e_file <- sprintf("Dependencies/Ebinds/Configuration_Files/Z_%003.f.csv", Z_choose)
unlink(e_file)
  
#Z_choose <- 24
# Grab neutral subshell binding energies
#assign K, L1 etc. notation
#grab all Eic's
#for(i in 1:100){
  Type <- "E2"
  Q <- 2000
  Z <- Z_choose
  bricc_string <- sprintf("./briccs -Z %d -g %.5f -L %s -a -w BrIccFO >> output.xml", Z, Q, Type)
  unlink("output.xml")
  system(bricc_string)
  ##Read data
  data <- xmlParse("output.xml")
  x <- xmlToList(data)
  len <- length(x) - 6

  ##Make the shells table
  Shells <- tibble(subshell = character(len), CCp = double(len), Eic = double(len) )
  for(k in 1:len){
    #k<-1
    # Shell <-  paste(unlist( x[6 + k]$PureCC$.attrs[[1]] ), collapse = "")
    # Shells$Shell[k] <- parse(text=Shell)
    #Shells$Shell[k] <- 2
    #k <- 2
    Shells$subshell[k]  <- x[5 + k]$PureCC$.attrs[[1]]
    Shells$Eic[k]  <- as.double(x[5 + k]$PureCC$.attrs[[2]])
    #Eic <- x[6 + k]$PureCC$.attrs[2]

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

## Assign subshell and compute Vi_q0 from Q-Eic
Shells$Vi_q0 <- Q-Shells$Eic
Shells <- filter(Shells, subshell %in% n_index$subshell)
Ebinds_n <- filter(Ebinds, Z==Z_choose)



## Assign the nh and nmax variables for the number of electrons in each subshell
for(i in 1:length(Shells$subshell)){
  Shells$nmax[i] <- n_index$nmax[which(n_index$subshell==Shells$subshell[i])[1]]
}

## Now smear out each observation of a subshell by the number of electrons in that subshell
ss_occ_vec <- rep(seq_len(nrow(Shells)), times = Shells$nmax)
ss_occ_vec <- ss_occ_vec[1:Z_choose]
Shells <- slice(Shells, ss_occ_vec)
Shells$Occ <- 1:Z_choose


i <- 1
while(i < length(Shells$subshell)){
  l <- length(which(n_index$subshell==Shells$subshell[i]))
  if(  (i+l) > length(Shells$subshell)  ){
    i <- length(Shells$subshell)
  } else {
  for(j in 1:l){
  # which(n_index$subshell==Shells$subshell[i]))[j]
    Shells$nh[i+j-1] <- n_index$nh[which(n_index$subshell==Shells$subshell[i])[j]]
  }
    i <- i + l   
  }
}

# Z_choose <- 23

#write.table(Shells, file=e_file, append=T, row.names=F, col.names=T,  sep=",")
atom <- Shells
Ebinds_n <- filter(Ebinds, Z==Z_choose)
atom$Z <- Z_choose

## Assign Ebinds and subshell colour and allot some missing values of Vi_q0
for(i in 1:nrow(Ebinds_n)){
  if( (is.na(Ebinds_n$IE[which(Ebinds_n$Occ == atom$Occ[i])]))  || is.null(Ebinds_n$IE[which(Ebinds_n$Occ == atom$Occ[i])]) ){
    atom$Ebind[i] <- NA
  } else {
    atom$Ebind[i] <- Ebinds_n$IE[which(Ebinds_n$Occ == atom$Occ[i])]
  }
  atom$shell_colours[i] <- Atomic_Colours$colour[which(Atomic_Colours$subshell==atom$subshell[i])]
  
  ## If the last Vi_q0 == NA or 0 then let it equal the first IE
  if(  (atom$Vi_q0[i] == 0) || (is.na(atom$Vi_q0[i]))   ){
    atom$Vi_q0[i] <- atom$Ebind[i]
  }
  
  ## Overall multiplier for ICC adjustment
  atom$multiplier[i] <- (atom$Ebind[i] / atom$Vi_q0[i]) * (atom$nh[i]/atom$nmax[i])
}

e_file <- sprintf("Dependencies/Ebinds/Configuration_Files/Z_%003.f.csv", Z_choose)
write.table(atom, file=e_file, append=T, row.names=F, col.names=T,  sep=",")

print(i)

}

ICC_Multipliers <- list.files(path = "Dependencies/Ebinds/Configuration_Files/",     # Identify all csv files in folder
                     pattern = "(Z_\\d+.csv)", full.names = TRUE) %>% 
  lapply(read_csv) %>%                                            # Store all files in list
  bind_rows                                             # Combine data sets into one data set 


for(i in 1:nrow(ICC_Multipliers)){
  ICC_Multipliers$Element[i] <- conv$Name[which(conv$Z == ICC_Multipliers$Z[i])]
}

##### Add in the bottom 10 elements
Ebinds_leq9 <- filter(Ebinds, Z<10) %>%
  select(Z, IE, Occ, Element) %>%
  rename(Ebind = IE)

ICC_Multipliers <- bind_rows(Ebinds_leq9, ICC_Multipliers)



write.table(ICC_Multipliers, file="Dependencies/ICC_Multipliers.csv", append=T, row.names=F, col.names=T,  sep=",")



