require(magrittr)
require(shiny)
require(shinyjs)
library(shinyWidgets)

runMetapopulation <- function(K = 20, N0 = 10, 
                              e = .1, c = .2, 
                              tmax = 1000,
                              extinction = FALSE,
                              colonization = FALSE){
  metapop.matrix <- matrix(ncol = K, nrow = tmax)
  metapop.matrix[1,] <- sample(c(rep(1, N0), rep(0, K-N0))) 
  
  for(i in 2:tmax){
    N.then <-  metapop.matrix[i-1,]
    f <- N.then / K
    if(extinction == FALSE)
      metapop.matrix[i,N.then == 1] <- rbinom(sum(N.then), size = 1, prob = 1-e) else 
        metapop.matrix[i,N.then == 1] <- rbinom(sum(N.then), size = 1, prob = (1-f)*e)
      if(colonization == FALSE)
        metapop.matrix[i,N.then == 0] <- rbinom(sum(N.then == 0), size = 1, prob = c) else 
          metapop.matrix[i,N.then == 0] <- rbinom(sum(N.then == 0), size = 1, prob = c*f)
  }
  
  return(list(mp = metapop.matrix, 
              Z =  runif(K) + 1i*runif(K),
              pars =  c(e = e, c = c),
              model = c(extinction, colonization)))
}

sim <- runMetapopulation(K = 100, extinction = FALSE, colonization = TRUE)

plotMP_timeseries <- function(sim){
  mp <- sim$mp
  e <- sim$pars["e"]
  c <- sim$pars["c"]
  
  rowSums(mp) %>% plot(type = "o", ylim = c(0,ncol(mp)), pch = 21, bg = "grey", 
                       xlab = "Time", ylab = "Occupied sites", col = "darkgrey")
  K <- ncol(mp)
  
  ex <- sim$model[1] == TRUE     # extinction process
  cz <- sim$model[2] == TRUE     # colonization process
  
  
  n.hat <- ifelse(!ex & !cz,  K*c/(c+e), 
                  ifelse(!ex & cz, max(0,K * (1-e/c)), 
                         ifelse(ex & !cz, min(K, K*c/e), NA)))
  
  linecols <- c("orange", "darkblue","darkred")
  abline(h = c(0, n.hat, ncol(mp)), col = linecols, lwd = 2)
  text(rep((nrow(mp)+1)/2, 3), c(0, n.hat, ncol(mp)), c("min", "theory", "max"), 
       col = linecols, pos = 3, adj = 1, xpd = NA, font = 2)
}

plotMP_timeseries(sim)


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
                 choices = list('Neutral' = FALSE, 'With rescue' = TRUE), 
                 selected = FALSE),
    
    radioButtons("colonization", "Colonization process:", 
                 choices = list('External' = FALSE, 'Internal' = TRUE),
                 selected = FALSE), 
    
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
    column(width = 7, style='padding:10px', plotOutput("plots")))
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
    
    par(mfrow = c(1,2), tck = 0.01, bty = "l", las= 1, 
        cex.lab = 1.5, mgp = c(1.25,.25,0))
    plotMP_timeseries(mysim)
    plot(mysim$Z, cex = 2, pch = 21, bg = mymp[nrow(mymp),], asp = 1, 
         xaxt = "n", yaxt = "n", xlab = "X", ylab = "Y")
  })
}


shinyApp(ui=ui, server=server)
