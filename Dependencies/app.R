# Sys.setenv("plotly_username" = "bw979")
# Sys.setenv("plotly_api_key" = "v7KDKKj5jCkkpKsid0kS")
# hello
library(rsconnect)
library(shiny)
library(tidyverse)
library(plotly)
library(DT)  ## for data tables
library(shinyjs)

source("Wins.R")
source("Rates_Plasma.R")
n_index <- read_csv("n_to_ICshell_Conversion.csv")
atomic_colours <- read_csv("Atomic_Colours.csv")

## Maxwell-Boltzman PD function
MB <- function(E, Te){
  F_E <- 2 * ((E/(pi))^(1/2)) * ((1/(Te))^(3/2)) * exp(-E / Te)
  return(F_E)
}

#facility <- filter(facility, Facility != c("XFELO", "SACLA", "EurXFEL", "LCLS")) %>%
facility <-  arrange(facility, Irradiance_Wcm2um2)

#### MAIN DATA TABLES
wap0 <- read_csv("All_RatesCalcd_Viva.csv")
Mo <- read_csv("Mo_Additions_Rates2.csv")
wap <- bind_rows(wap0, Mo)
wap <- select(wap, -Rneec_Pl, -Rneec_SEBIT) %>%
  filter(Z>9)
wap_N0 <- filter(wap, Ei==0)
wap_Ni <- filter(wap, Ei!=0)

# isomer$N <- isomer$M - isomer$Z
# isomer$Th_numeric <- pblapply(isomer$Thalf, as_Numeric_Th_iso)
# isomer$Th_string <- pblapply(isomer$Thalf, as_Type_Th_iso)
#source("Wins_initialize.R")
####################################################################################################
### UI #############################################################################################
####################################################################################################
ui <- fluidPage(
  useShinyjs(),
  navbarPage("Select:",
             ## tabPanel("name", content)
             tabPanel("Find Nuclide",
                      
                      # sidebarPanel(
                      #   sliderInput(inputId = "energy0",label = "Max Impact Energy (keV)", value = 5000, min = 1, max = 5000)
                      # ),   
                      numericInput("energy0", "Max Nuclear Transition Energy / keV", value = 4000, min = 1, max = 4000),  
                      #mainPanel(   
                      ## PLot output
                      #plotOutput(outputId = "MasterPlot"), 
                      ## Data Table output
                      DT::dataTableOutput("MasterTable"),  
                      DT::dataTableOutput("FACILITY")
                      
                      #tableOutput("MasterTable")  
                      #)   
             ),
             tabPanel("Explore Nuclide",
                     #always need an input from the user
                     tags$head(tags$style(".shiny-notification {position: fixed; top: 10% ;left: 50%}")),
                      sidebarPanel(
                        p("Note: can take several minutes to calculate NEEC rates and produce plots. There may be errors displayed whilst the NEEC rates are being calculated and plotted"),
                        #p("ne, Te_max and repetition rate are defined by the facility"),
                        
                        # Choice of input type
                        radioButtons(inputId="In_Type", label="Input type?", 
                                     choices=c("Astrophysical","Facility"), selected = "Facility"),
                        #selectInput("In_Type", "Input Type:", choices=c("Astro", "Facility"), selected="Astro"),
                        
                        # submit button
                        actionButton("Go", label = "Go"),#
                        
                        
                        
                        # sliderInput(inputId = "energy",label = "Max Impact Energy (keV)", value = 5000, min = 1, max = 5000),
                        textInput("AX", "Which nuclide (element in caps)", "57FE"),
                        selectInput("variable", "Initial Nuclear Level:", choices=list("Ground State" = "gs", "Isomeric" = "iso"), selected="Ground State"),
                        uiOutput("Nuc_Transitions"),
                        selectInput("Out_Type", "Output Type:", choices=c("Number", "Rate"), selected="Number"),
                        #actionButton("goButton", "Go!", class = "btn-success"),
                        uiOutput("Te_Input"),
                        uiOutput("ne_Input"),
                        
                        uiOutput("Ep_J"),  
                        uiOutput("Rf_um"),
                        uiOutput("tl_fs"),
                        uiOutput("lambda_nm"),
                        
                        
                        #### MACROSCOPIC
                        #selectInput("facil", "Laser Facility:", choices=facility$Facility, selected="ELI-NP"),
                        uiOutput("facil_ui"),
                        uiOutput("Print_ne_max"), 
                        uiOutput("Print_Te_max"), 
                        uiOutput("Print_RepRate"),
                        uiOutput("Te_facility"),
                      
                        #sliderInput(inputId = "Te_fine",label = "Electron Temperature (Fine Tune) (keV)", value = 5, min = 0, max = 10),
                        sliderInput(inputId = "tplasma",label = "Plasma Lifetime  (ps)", value = 100, min = 0.1, max = 1000),
                        #numericInput("RepRate", "Laser Repetition Rate (Hz):", value = 0.0167, min = 0.0000001, max = 100000000),
                        sliderInput(inputId = "texp",label = "Length of Experiment (Weeks):", value = 1, min = 0.1, max = 5),
                        #plotlyOutput("side_chart"),
                        # checkboxGroupInput("variable", "Initial State of Target:",
                        #                   c("Isomeric" = "iso",
                        #                     "Ground_State" = "gs"),
                        #                    selected = "gs"
                        #)
                        
                       

                        
                        #DT::dataTableOutput("IsoInfo") 
                        plotlyOutput("side_chart"),
                        
                        # ### Optimiser Parameters
                        textInput("E1", "Input Lower Electron Temperature (keV)"),
                        textInput("E2", "Input Upper Electron Temperature (keV)"),
                        textInput("Estep", "Stepsize"),
                        # submit button
                        actionButton("Optimise", label = "OPTIMISE TEMPERATURE")#
                        
                        
                        #checkboxInput("iso", "Isomeric", FALSE)
                        # verbatimTextOutput("value")
                      ),
                      #always need an output type to pass to server
                      mainPanel(
                        
                        plotlyOutput("bar_chart"),
                        #p("Total Number of NEECs in this scenario:"),
                        plotlyOutput("bar_chart2"),
                        plotlyOutput("CSD_Chart"),
                        plotlyOutput("Optim_Chart"),
                        DT::dataTableOutput("FilteredTable")  
                      )  
             )
  )
)
WapTable <- select(wap,-Ji_double, -Jf_double, -Thi_numeric, -Thf_numeric, -M, -El, -N, -Egam, -FL) 
#HANDY
#filter(wap_Ni, !is.na(S), S > 0)







####################################################################################################
## SERVER ##########################################################################################
####################################################################################################
server <- function(input, output) {
  
  ##### SELECT FACILITY
  output$facil_ui <- renderUI({
    selectInput("facil", "Laser Facility:", choices=facility$Facility, selected="ELI-NP")
  })
  
  # ### Optimiser Parameters
  ## REACTIVE input for Te
  E1_inp <- reactive({
    as.double(input$E1) #+ input$Te_fine
  })
  
  E2_inp <- reactive({
    as.double(input$E2) #+ input$Te_fine
  })
  
  Estep_inp <- reactive({
    as.double(input$Estep) #+ input$Te_fine
  })
 

  
  Ee_max <- 5000
  
  ### Rename variables
  output$MasterTable <- DT::renderDataTable({
    WapTable2 <- rename(WapTable, "T1/2 initial" = "Thi", "T1/2 final" = "Thf","Ei keV" = "Ei", "Ef keV" = "Ef",  "Capture Level Configuration" = "At_GS_Config", "Charge State" = "CS", "Ee keV" = "Ee", "Q keV" = "Q", "Atomic Ebind keV" = "Ebind")
    DT::datatable(WapTable2)
  })
  
  ### List the possible nuclear transitions
  output$Nuc_Transitions <- renderUI({
    wap_filtered <- wap_f()
    Trans <- as.character(unique(wap_filtered$Q))
    selectInput("Transitions", "Select Nuclear Transition", c("All", Trans))
    
  })
  
  ## INPUT Te directly
  output$Te_Input <- renderUI({
    #sliderInput(inputId = "Te",label = Te_string_max, value = Te_Max, min = 1, max = Te_Max)
    textInput(inputId = "Te", label = "Enter Electron Temperature (keV)", value = 1000, width = NULL, placeholder = NULL)
  })
  
  
  ## INPUT ne directly
  output$ne_Input <- renderUI({
    #sliderInput(inputId = "Te",label = Te_string_max, value = Te_Max, min = 1, max = Te_Max)
    textInput(inputId = "ne", label = "Enter Plasma Number density (e-'s cm^-3)", value = 1E20, width = NULL, placeholder = NULL)
  })
  
  
#### OBSERVE EVENT: Change in input type radio button  
observeEvent(input$In_Type, { 
  
###### ASTROPHYSICAL ######
if(input$In_Type == "Astrophysical") {

  #hide facility selection and plasma sliders
  hide(id = "facil")
  hide(id = "tplasma")
  hide(id = "texp")
  
  show("Te")
  show("ne")
  
  #### Print maximal facility Te value
  output$Print_Te_max <- NULL
  
  #### Print Maximal facility ne value
  output$Print_ne_max <- NULL
  
  #### Print Maximal Rep Rate
  output$Print_RepRate <- NULL
   
  ## uiOutput("Ep_J") 
  output$Ep_J <- NULL
  #hide("Epulse_J")
  
  ## uiOutput("Rf_um")
  output$Rf_um <- NULL
  #hide(id="Rfocal_um")
  
  #uiOutput("tl_fs")
  output$tl_fs <- NULL
  #hide(id="tlaser_fs")
  
  #uiOutput("lambda_nm")
  output$lambda_nm <- NULL
  
 }
###### END ASTROPHYSICAL ######
  
  
###### FACILITY ######
if(input$In_Type == "Facility") {
  
  #hide facility selection and plasma sliders
  show(id = "facil")
  show(id = "tplasma")
  show(id = "texp")
  
  hide(id = "Te")
  hide(id = "ne")
  
  #### Print maximal facility Te value
  output$Print_Te_max <- renderUI({
    Te_facil_max <- signif(filter(facility, Facility==input$facil)$Te_chosen_keV, 3)
    Te_string_max <- paste("Maximal electron temperature at", input$facil, "is", Te_facil_max, "keV")
    p(paste("Maximal electron temperature at", input$facil, "is", Te_facil_max, "keV"))
  })
  
  #### Print Maximal facility ne value
  output$Print_ne_max <- renderUI({
    ne_facil_max <- signif(filter(facility, Facility==input$facil)$ne, 3)
    p(paste("Maximal number density at", input$facil, "is", ne_facil_max, "e-'s cm^-3"))
  })
  
  
  #### Print Maximal Rep Rate
  output$Print_RepRate <- renderUI({
    RepRate <- signif(filter(facility, Facility==input$facil)$Max_Repetition_Rate_Hz, 3)
    p(paste("Maximal laser rep rate at", input$facil, "is", RepRate, "Hz"))
  })
  
  ## uiOutput("Ep_J") 
  output$Ep_J <- renderUI({
    Ep <- signif(filter(facility, Facility==input$facil)$Pulse_Energy_J, 3)
    l <- p(paste("Enter Pulse Energy, facility maximum is ", Ep, " J"))
    textInput(inputId = "Epulse_J", label = "Enter Pulse Energy (J)", value = Ep, width = NULL, placeholder = NULL)
  })
  #Ep_inp <- reactive({ as.double(input$Epulse_J) })
  
  ## uiOutput("Rf_um")
  output$Rf_um <- renderUI({
    Rf <- 0.5 * signif(filter(facility, Facility==input$facil)$Spot_Diameter_um, 3)
    l <- p(paste("Enter Spot Size, facility minimum is ", Rf, " um"))
    textInput(inputId = "Rfocal_um", label = l, value = Rf, width = NULL, placeholder = NULL)
  })
  #Rf_inp <- reactive({ as.double(input$Rfocal_um) })
  
  #uiOutput("tl_fs")
  output$tl_fs <- renderUI({
    tl <- signif(filter(facility, Facility==input$facil)$Pulse_Duration_fs, 3)
    l <- p(paste("Enter Pulse duration, facility minimum is ", tl, " fs"))
    textInput(inputId = "tlaser_fs", label = l, value = tl, width = NULL, placeholder = NULL)
  })
  #tl_inp <- reactive({ as.double(input$tlaser_fs) })
  
  #uiOutput("lambda_nm")
  output$lambda_nm <- renderUI({
    lambda <- signif(filter(facility, Facility==input$facil)$Wavelength_nm, 3)
    #l <- p(paste("Enter Pulse duration, facility minimum is ", tl, " fs"))
    textInput(inputId = "wavelength_nm", label = "Enter Laser Wavelength (nm)", value = lambda, width = NULL, placeholder = NULL)
  })
  #lambda_inp <- reactive({ as.double(input$wavelength_nm) })
  # wf2$fq[i] <- compute_Plasma_Rate(wf2[i,], TRUE, facility[j,], (input$tplasma)*1E-12, "params", 1,1, input$Epulse_J ,input$Rfocal_um , input$tlaser_fs, input$wavelength_nm, "fq", "p") 
  # wf2$Rate_Pl[i] <- compute_Plasma_Rate(wf2[i,], TRUE, facility[j,], (input$tplasma)*1E-12, "params", 1,1, input$Epulse_J ,input$Rfocal_um , input$tlaser_fs ,input$wavelength_nm, "Rate", "p") 
  
  
}  
  ##### END FACILITY ######
}) ###### END OBSERVE RADIO BUTTON EVENT #####
 
## REACTIVES
Te_inp <- reactive({ as.double(input$Te)})
ne_inp <- reactive({ as.double(input$ne) })

Ep_inp <- reactive({ as.double(input$Epulse_J) })
Rf_inp <- reactive({ as.double(input$Rfocal_um) })
tl_inp <- reactive({ as.double(input$tlaser_fs) })  
lambda_inp <- reactive({ as.double(input$wavelength_nm) })


## Reactive for isomer or ground initial state
wap_f <- reactive({
  if(input$variable == "iso"){
    wf <<- filter(wap_Ni, Ee < Ee_max, AX == input$AX)
  }
  else if(input$variable == "gs"){
    wf <<- filter(wap_N0, Ee < Ee_max, AX == input$AX)
  }
  return(wf)
}) 


##### Reactive to recalculate NEEC rates HERE
wap_f2 <- reactive({
  ## filter by chosen nuclear transition or 'all' first
  if(input$Transitions == "All"){  
    wf2 <- wap_f()
  } else {
    wf2 <- filter(wap_f(), Q == input$Transitions)
  }
  ## Set MB distribution values
  if(is.null(wf2$Ee)) {print("NO ISOMERS")} else {
    Evals <- c(seq(0, max(wf2$Ee), 0.1 ))
    MBvals <- mapply(MB, Evals, Te_inp()) 
  }  
  
  if(input$In_Type == "Facility"){
    ## GET FACILITY INDEX
    j <- which(facility$Facility == input$facil)
    
    for( i in 1:length(wf2$AX) ) {
      withProgress(message = "Calculating Rates and Reading CSD's", value = (i/length(wf2$AX)), { 
        #wap_filtered$Unity_Plasma_Rate[i] <- compute_Plasma_Rate(wap_filtered[i,], facility[j,], "Te", Te_inp(), "Rate", "p")
        #compute_Plasma_Rate(Wap_filtered, facility[j,], "Te", "Rate", "p")
        #wf2$S[i] <- wf2$S[i] * wf2$Gamma_Scale_Factor[i]
        ####compute_Plasma_Rate <- function(Wins_all, FACIL_bool, facility_in, tau, INPUT_Te_Type, Te_input, ne_input, Ep_J, Rfoc_um, tlas_fs, Wavelen_nm, OUTPUT, OUTPUT2)
        wf2$fq[i] <- compute_Plasma_Rate(wf2[i,], FALSE, facility[j,], (input$tplasma)*1E-12, "Params", 1000,1E20, Ep_inp() , Rf_inp() , tl_inp(), lambda_inp(), "fq", "p") 
        wf2$Rate_Pl[i] <- compute_Plasma_Rate(wf2[i,], FALSE, facility[j,], (input$tplasma)*1E-12, "Params", 1000, 1E20, Ep_inp() , Rf_inp() , tl_inp(), lambda_inp(), "Rate", "p") 
        wf2$Npw_Pl_S[i] <- facility$Max_Repetition_Rate_Hz[j]*input$tplasma*(1E-12)*input$texp*60*60*24*7 * wf2$Rate_Pl[i]
      })
    }
  } else if(input$In_Type == "Astrophysical"){
    ## GET FACILITY INDEX
    j <- 1
    for( i in 1:length(wf2$AX) ) {
      withProgress(message = "Calculating Rates and Reading CSD's", value = (i/length(wf2$AX)), {
        #wap_filtered$Unity_Plasma_Rate[i] <- compute_Plasma_Rate(wap_filtered[i,], facility[j,], "Te", Te_inp(), "Rate", "p")
        #compute_Plasma_Rate(Wap_filtered, facility[j,], "Te", "Rate", "p")
        #wf2$S[i] <- wf2$S[i] * wf2$Gamma_Scale_Factor[i]
        wf2$fq[i] <- compute_Plasma_Rate(wf2[i,], FALSE, facility[j,], (input$tplasma)*1E-12, "Te", 1000*Te_inp(), ne_inp(), 1,1,1,1, "fq", "p")
        wf2$Rate_Pl[i] <- compute_Plasma_Rate(wf2[i,], FALSE, facility[j,], (input$tplasma)*1E-12, "Te", 1000*Te_inp(), ne_inp(), 1,1,1,1, "Rate", "p")
        wf2$Npw_Pl_S[i] <- facility$Max_Repetition_Rate_Hz[j]*input$tplasma*(1E-12)*input$texp*60*60*24*7 * wf2$Rate_Pl[i]
      })
    }
  }
  return(wf2)  
})







##### Do the plots and recalculate when Go button is pressed  
observeEvent(input$Go, {
  


 
  
  
  
  ###### ALL PLOTS BELOW HERE
  output$bar_chart <- renderPlotly({
    input$goButton
    wap_filtered2 <- wap_f2()
    
    if(is.null(wap_filtered2$Ee)) {print("NO ISOMERS")} else {
      Evals <- c(seq(0, max(wap_filtered2$Ee), 0.1 ))
      MBvals <- mapply(MB, Evals, Te_inp()) 
    }  
   
    
    #### RESONANCE STRENGTH CHART
    plot_ly(wap_filtered2, x = ~Ee, y = ~S, type = 'bar', showlegend=TRUE, name=~shell) %>% #, marker = list(color = label_colours)) %>%
      layout(yaxis = list(type = "log", exponentformat = "E", title = "NEEC Resonance Strength / beV"), #range = c(1E-10, max(wap_filtered$S))),
             xaxis = list(tickangle = -45, title="Electron Impact Energy / keV", range = c(min(wap_filtered2$Ee), max(wap_filtered2$Ee))),
             legend=list(name=~shell, showlegend=TRUE)
             # yaxis2 = list(
             #   tickfont = list(color = "red"),
             #   overlaying = "y",
             #   side = "right",
             #   title = "MB" ,
             #   showlegend = F)
        
      )  %>%
      add_lines(x = ~Evals, y = ~MBvals, type = 'scatter', mode = 'lines', color="#DE3163", showlegend=T, name="F(E)") 
   
    
    
    #plot_ly(  isomer, x = ~1:length(isomer$M), y = ~input$isomer_energy, type = 'bar')  
    #hist(rnorm(input$isomer_energy))
  })
  
  #### YIELD CHART
  output$bar_chart2 <- renderPlotly({
  
    wap_filtered2 <- wap_f2()
    if(input$Out_Type == "Number"){
    wap_filtered2 <- filter(wap_filtered2, !is.na(Npw_Pl_S))
    N_count <- signif(sum(wap_filtered2$Npw_Pl_S), 5)
    pr <- paste("Sum Yield:", N_count)
    
    plot_ly(wap_filtered2, x = ~Ee, y = ~Npw_Pl_S, type = 'bar', name =~shell) %>% layout( # , marker=list(color=wap_filtered$shell_colour))
      yaxis = list(type = "log", exponentformat = "E", title = "Number of NEEC's"), #range = c(1E-10, max(wap_filtered$Npw_Pl_S))),
      xaxis = list( title= "Electron Impact Energy / keV", range = c(min(wap_filtered2$Ee), max(wap_filtered2$Ee))),
      showlegend = TRUE,
      annotations = list(text = pr, showarrow=FALSE )
      )
    } else {
      wap_filtered2 <- filter(wap_filtered2, !is.na(Rate_Pl))
      N_count <- signif(sum(wap_filtered2$Rate_Pl), 5)
      pr <- paste("Sum RATE (s^-1):", N_count)
      
      plot_ly(wap_filtered2, x = ~Ee, y = ~Rate_Pl, type = 'bar', name =~shell) %>% layout( # , marker=list(color=wap_filtered$shell_colour))
        yaxis = list(type = "log", exponentformat = "E", title = "NEEC Rate / s^-1"), #range = c(1E-10, max(wap_filtered$Npw_Pl_S))),
        xaxis = list( title= "Electron Impact Energy / keV", range = c(min(wap_filtered2$Ee), max(wap_filtered2$Ee))),
        showlegend = TRUE,
        annotations = list(text = pr, showarrow=FALSE )
      )
    }
    
  })
  
  output$CSD_Chart <- renderPlotly({
   
    wap_filtered2 <- wap_f2()
    
    plot_ly(wap_filtered2, x = ~Occ, y = ~fq, type="scatter", mode = 'markers+lines', color="yellow") %>%
      layout(
        yaxis=list(title="Pq"),
        xaxis=list(title="Number of Electrons on Ion - Pre Capture")
      )
    
  })
  
  ### SIDE CHART
  output$side_chart <- renderPlotly({
    if(input$variable == "gs" && input$variable == "iso"){
      wap_filtered <- bind_rows(filter(wap_N0, Ee < Ee_max, AX == input$AX) , filter(wap_Ni, Ee < Ee_max, AX == input$AX) ) 
    }
    else if(input$variable == "iso"){
      wap_filtered <- filter(wap_Ni, Ee < Ee_max, AX == input$AX)  
    }
    else if(input$variable == "gs"){
      wap_filtered <- filter(wap_N0, Ee < Ee_max, AX == input$AX)  
    }
    
    # Bar graph showing isomer energies
    plot_ly(wap_filtered, x =~AX, y =~Ei, type="bar", opacity=1.0, size=4, name=~Thi, title="Initial Energies", text = ~Thi) %>%
      layout(
        yaxis=list(exponentformat='E', title="Isomer Energy / keV")
      )
  })
  
  
 
  
}) ### end of  observe event GO brackets  
  

observeEvent(input$Optimise, {
  #### Optim_Chart
  output$Optim_Chart <- renderPlotly({
    #wap_filtered2 <- wap_f2()
    
    #xvals ... input e's
    xvals <- c(seq(E1_inp(), E2_inp(), Estep_inp()))
    
    #yvals
    ##yvalue array:
    N_count <- c(rep(0, length(xvals)))
    ## filter by chosen nuclear transition or 'all' first
    if(input$Transitions == "All"){  
      wap_filtered2 <- wap_f()
    } else {
      wap_filtered2 <- filter(wap_f(), Q == input$Transitions)
    }
    
    j <- 1
    
    for(k in 1:length(xvals)) {
      # ## Set MB distribution values
      # if(is.null(wf2$Ee)) {print("NO ISOMERS")} else {
      #   #Evals <- c(seq(0, max(wf2$Ee), 0.1 ))
      #   MBvals <- mapply(MB, xvals[k], Te_inp()) 
      # }  
      
      for( i in 1:length(wap_filtered2$AX) ) {
        withProgress(message = "Calculating Rates and Reading CSD's", value = (i/length(wap_filtered2$AX)), {
          #wap_filtered$Unity_Plasma_Rate[i] <- compute_Plasma_Rate(wap_filtered[i,], facility[j,], "Te", Te_inp(), "Rate", "p")
          #compute_Plasma_Rate(Wap_filtered, facility[j,], "Te", "Rate", "p")
          #wf2$S[i] <- wf2$S[i] * wf2$Gamma_Scale_Factor[i]
          #wf2$fq[i] <- compute_Plasma_Rate(wf2[i,], FALSE, facility[j,], (input$tplasma)*1E-12, "Te", 1000*Te_inp(), ne_inp(), 1,1,1,1, "fq", "p")
          wap_filtered2$Rate_Pl[i] <- compute_Plasma_Rate(wap_filtered2[i,], FALSE, facility[j,], (input$tplasma)*1E-12, "Te", 1000*xvals[k], ne_inp(), 1,1,1,1, "Rate", "p")
          #wf2$Npw_Pl_S[i] <- facility$Max_Repetition_Rate_Hz[j]*input$tplasma*(1E-12)*input$texp*60*60*24*7 * wf2$Rate_Pl[i]
        })
      }
      wap_filtered2 <- filter(wap_filtered2, !is.na(Rate_Pl))
      ## Append this to the yvalue array
      N_count[k] <- signif(sum(wap_filtered2$Rate_Pl), 5)
    }
    #return(wf2)  
    
    
    
    
    
    ##Energy Range
    #E_range <- input$E2 - input$E1
    
    
    #wap_filtered2 <- filter(wap_filtered2, !is.na(Rate_Pl))
    #N_count <- signif(sum(wap_filtered2$Rate_Pl), 5)
    
    plot_ly(x = xvals, y = N_count, type="scatter", mode = 'markers+lines', color="blue") %>%
      layout(
        yaxis=list(title="Total NEEC Rate (s^-1)"),
        xaxis=list(title="Temperature (keV)")
      )
    
  })
  
})
  
  
  ##### DATA TABLES ########
  ### DATA TABLE 1
  output$FilteredTable <- DT::renderDataTable({
    DT::datatable(wap_f2())
  })
  
  ### DATA TABLE 2
  output$FACILITY <- DT::renderDataTable({
    DT::datatable(facility)
  })
  
  
}

shinyApp(ui = ui, server = server)

