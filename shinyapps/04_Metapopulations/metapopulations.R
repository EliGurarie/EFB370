require(magrittr)
require(shiny)
require(shinyjs)
library(shinyWidgets)

source("metapopulation_functions.R")

eval <- FALSE
if(eval){

  e <- .3; c <- .4 
  sim <- runMetapopulation(K = 100, 
                           extinction = "constant", 
                           colonization = "external",
                           N0 = 10, e = e, c = c, 
                           tmax = 1000)
  plotMP_timeseries(sim)
  abline(h = (1-e/c)*100)
  abline(h = (c/e)*100)
}

ui <- fluidPage(
  
  useShinyjs(),
  setBackgroundColor("MidnightBlue"),
  tags$style(
    "*, div {
      font-family: Consolas;
    }"
  ),
  
  h1("Metapopulation Modeller 3000XC", style = "color: HotPink"),
  
  sidebarPanel(
    
    h3("Process model"),
    
    radioButtons("extinction", label = "Extinction process:",
                 choices = list('Neutral' = "constant", 'With rescue' = "rescue"), 
                 selected = "constant"),
    
    radioButtons("colonization", "Colonization process:", 
                 choices = list('External' = "external", 'Internal' = "internal"),
                 selected = "external"), 
    
    h3("Metapopulation setup"),
    
    numericInput("K", "Total number of sites (K)", 
                 value = 20, min = 1, step = 1),
    
    numericInput("N0", "Initial number of local populations (N0)", 
                 value = 5, min = 1, step = 1),
    
    sliderInput(inputId = "e",
                label = "Extinction probability",
                min = 0, max = 1, value = 0.1),
    
    sliderInput(inputId = "c",
                label = "Colonization probability",
                min = 0, max = 1, value = 0.1),
    
    numericInput(inputId = "tmax", 
                 label = "Duration:", 
                 value = 2, min = 2, max = 400),
    actionButton("go", "Refresh simulation")),
  
  fluidRow( 
    column(width = 7, height = 7, style='padding:10px', plotOutput("plots", width = "600px", height = "500px")))
)

server <- function(input, output) {
  sim <- eventReactive(input$go,
                       runMetapopulation(K = input$K, N0 = input$N0,
                                         e = input$e, c = input$c, tmax = 400, 
                                         extinction = input$extinction,
                                         colonization = input$colonization))
  output$plots <- renderPlot({
    mysim <- sim()
    mymp <- mysim$mp[1:input$tmax,]
    mysim$mp <- mymp
    
    par(tck = 0.01, bty = "l", las= 1, 
        cex.lab = 1.5, mgp = c(1.25,.25,0))
    
    layout(rbind(c(1,1,2,2),c(4,3,3,5)))
    plotMP_timeseries(mysim)
    plot(mysim$Z, cex = 2, pch = 21, bg = mymp[nrow(mymp),], asp = 1, 
         xaxt = "n", yaxt = "n", xlab = "X", ylab = "Y")
    plotMP_theory(mysim)
  })
}


shinyApp(ui=ui, server=server)
