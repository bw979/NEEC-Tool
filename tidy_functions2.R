# tidies the "J" column in any NNDC data.frame, outputs filtered df with options for tentative assignments
tidy_Jpi <- function(input){
  input <- input %>% mutate(Tentative_Jpi = logical(1))
  input <- input %>% mutate(Multiple_Jpi = logical(1))
  input <- input %>% filter(!is.na(J))
  #replace with clean
  for(i in 1:length(input$J)) {
    print(input$AX[i])
    if (str_detect(input$J[i], "3/2 5/2") || str_detect(input$J[i], "3/2\n5/2")) {
      input$J[i] <- "3/2"
    }
    if(!is.na(input$J[i]) &&  str_detect(input$J[i], "[\\[()\\]]")   ) {
      input$Tentative_Jpi[i] = TRUE
      input$J[i] <- str_remove_all(input$J[i], "[\\[()\\]]")
    }
    if(!is.na(input$J[i]) && str_detect(input$J[i], "[A-z\\,\\&\\:]") ) {
      input$Multiple_Jpi[i] = TRUE
    }
    if(!is.na(input$J[i]) && input$J[i] == "") {input$J[i] = NA} 
    
    ## Get rid of the "+"
    if(!is.na(input$J[i]) && nchar(input$J[i])==1 &&  str_detect(input$J[i], "\\+")   ) {
      input$Multiple_Jpi[i] = TRUE
      #input$J[i] <- str_remove_all(input$J[i], "[\\[()\\]]")
    }
  }
  #modified input up to here
  output <- input %>% #filter(Multiple_Jpi == FALSE) %>%
   filter(!is.na(J))
  return(output)
}
