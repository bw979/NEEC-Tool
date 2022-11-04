###############################
######### CALCULATE ICC #######
print("Calculate ICC using BRICCS")

#####ICC's using map function
AnNi$ICC_tot <- pmap(list(AnNi$Type, AnNi$Z, AnNi$Q, "tot"), Compute_ICC_singles )
AnNi$ICC <- pmap(list(AnNi$Type, AnNi$Z, AnNi$Q, "p"), Compute_ICC_singles )
#AnNi$Ebind_q0 <- pmap(list(AnNi$Type, AnNi$Z, AnNi$Q, "Eic"), Compute_ICC_singles )

unlink(Name_AnN0_Rates)
write.table(AnN0, file=Name_AnN0_Rates, append=T, row.names=F, col.names=T,  sep=",")
AnN0 <- read_csv(Name_AnN0_Rates)

### Filter only the B vals that are known
# length(unique(AnN0$AX))
# AnN0_All_B <- AnN0
#AnN0 <- filter(AnN0, !is.na(B), !is.na(IC), Type != "E0")
#length(unique(An2$AX))

##### Have to extract J value and then evaulate in a new numeric column
AnN0$Jf_double <- double(1)
AnN0$Ji_double <- double(1)

print("Compute Ar and S")

for(i in 1:length(AnN0$Ar)){
  #i<-1
  print("Compute Ar and S")
  print(i)
  if(AnN0$Ji[i]=="+") AnN0$Ji_double[i] = NA
  if(AnN0$Jf[i]=="+") AnN0$Jf_double[i] = NA
  if(AnN0$Ji[i]=="-") AnN0$Ji_double[i] = NA
  if(AnN0$Jf[i]=="-") AnN0$Jf_double[i] = NA
  
  if(AnN0$Ji[i]!="+" && AnN0$Jf[i]!="+" && AnN0$Ji[i]!="-" && AnN0$Jf[i]!="-"){
    
    # Make the Jpi strings into just J doubles
    # AnN0$Jf_double[i] <- paste(unlist(str_extract_all(AnN0$Jf[i], "[0-9/]")), collapse = "")
    # AnN0$Jf_double[i] <- eval(parse(text=AnN0$Jf_double[i]))
    
    #AnN0$Jf_double[i] <- paste(unlist(str_extract_all(AnN0$Jf[i], "[0-9/]")), collapse = "")
    AnN0$Jf_double[i] <- eval(parse(text=paste(unlist(str_extract_all(AnN0$Jf[i], "[0-9/]")), collapse = "")))
    
    # AnN0$Ji_double[i] <- paste(unlist(str_extract_all(AnN0$Ji[i], "[0-9/]")), collapse = "")
    # AnN0$Ji_double[i] <- eval(parse(text=AnN0$Ji_double[i]))
    AnN0$Ji_double[i] <- eval(parse(text=  paste(unlist(str_extract_all(AnN0$Ji[i], "[0-9/]")), collapse = "")  ))
    
    # AnN0$Jf_double[i] <- as.double(as.character(AnN0$Jf_double[i]))
    # AnN0$Ji_double[i] <- as.double(as.character(AnN0$Ji_double[i]))
  }
  if(  !is.na(AnN0$B[i]) && !is.na(AnN0$Ji_double[i]) && !is.na(AnN0$Jf_double[i]) && !is.na(AnN0$ICC[i]) && AnN0$Type[i] != "E0" && AnN0$Type[i] != "Large_L") {
    ######Compute_Ar
    AnN0$Ar[i] <- compute_Ar(AnN0$Type[i], AnN0$M[i], AnN0$Q[i], AnN0$B[i])
    #Compute_S Jd, Jat,Ji,ICC, Ar, Etrans_keV, Ebind_keV
    AnN0$S[i]  <-  compute_S(AnN0$Jf_double[i], AnN0$Jat[i], AnN0$Ji_double[i], AnN0$ICC[i], AnN0$Ar[i], AnN0$Q[i], AnN0$Ebind[i])
  }
}