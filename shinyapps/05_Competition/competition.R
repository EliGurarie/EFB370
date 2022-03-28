require(magrittr)
require(shiny)
require(shinyjs)
library(shinyWidgets)

runCompetition <- function(r.n, r.m, K.n, K.m, N0, M0, alpha, beta, tmax){
  
  N <- matrix(0, ncol = K.n, nrow = tmax)
  M <- matrix(0, ncol = K.m, nrow = tmax)
  
  N[1,sample(1:K.n, N0)] <- 1
  M[1,sample(1:K.m, M0)] <- 1
  
  for(i in 2:tmax){
    
    # N and M are vectors of occupancy
    # n and m are counts
    
    N.occupied <- which(N[i-1,] == 1)
    n.new <- sum(N[i-1,])*r.n
    N.potential <- unique(c(N.occupied, sample(1:K.n, n.new)))
    
    M.occupied <- which(M[i-1,] == 1)
    m.new <- sum(M[i-1,])*r.m
    M.potential <- unique(c(M.occupied, sample(1:K.m, m.new)))
    
    n.removed <- rpois(1, alpha * length(M.occupied) * length(N.potential) /  K.n)
    m.removed <- rpois(1,  beta * length(N.occupied) * length(M.potential) /  K.m)
    
    N[i,sample(N.potential, max(length(N.potential) - n.removed, 0))] <- 1
    M[i,sample(M.potential, max(length(M.potential) - m.removed, 0))] <- 1
  }
  
  
  Z.m <- runif(K.m, 0, 1) + 1i * runif(K.m, 0, 1)
  Z.n <- runif(K.n, 0, 1) + 1i * runif(K.n, 0, 1)
  
  return(list(M = M, N = N, Z.m = Z.m, Z.n = Z.n))
}

plotMN <- function(t, sim){
  par(mar = c(3,3,2,2), bty = "l", mgp = c(1.5,.25,0),
      tck = 0.01, cex.lab = 1.5, cex.axis = 1.25)
  
  N <- sim$N
  M <- sim$M
  Z.n <- sim$Z.n
  Z.m <- sim$Z.m
  
  N.time <- apply(N, 1, sum)
  M.time <- apply(M, 1, sum)
  
  layout(rbind(c(1,1), 2:3))
  
  bgs <- c("pink", "lightblue")
  cols <- c("darkorange", "darkblue")
    
  plot(c(0,min(2*t, nrow(N))), c(0, max(N.time, M.time)), type = "n", ylab = "population", xlab = "time")
  points(N.time[1:t], type = "o", pch = 21, bg = "pink", col = cols[1], cex = 2)
  points(M.time[1:t], type = "o", pch = 21, bg = "lightblue", col = cols[2], cex = 2)

  legend("right", pch = 21, pt.bg = bgs, col = cols, legend = c("Squirlicorn", "Pegamunk"), cex = 2)
  
  plot(0,0,xlim = c(0,1), ylim = c(0,1), asp = 1, type = "n", 
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "o")
  points(Z.m[M[t,] == 1], bg = bgs[1], col = cols[1], pch = 21, cex = 2)
  points(Z.n[N[t,] == 1], bg = bgs[2], col = cols[2], pch = 21, cex = 2)

  plot(N.time[1:t], M.time[1:t], xlim = c(0,K.n), ylim = c(0,K.m), 
       xlab = "Squirlicorn", ylab = "Pegamunk", col = "purple", lwd = 1.5, bg = "green", 
       type = "o", pch = 21, cex = 2, asp = 1)
  
  abline(v = K.m, col = "grey", lwd = 2, lty = 3)
  abline(h = K.n, col = "grey", lwd = 2, lty = 3)
  
   #abline(K.n/alpha, -1/alpha, col = "darkorange", lwd = 2)
   #abline(K.m, -beta, col = "darkblue", lwd = 2)
}



ui <- fluidPage(
  
  useShinyjs(),
  setBackgroundColor("Black"),
  tags$style(
    "*, div {
      font-family: 'Comic Sans MS', 'Comic Sans', cursive;
    }"
  ),
  
  h1("Squirlicorn vs. Pegamunk: The Final Chapter", style = "color: HotPink"),
  
  
  sidebarPanel(
    img(src='squirlicornvspegamunk_small.png', align = "center"),
               
    h3("Squirlicorn parameters:"),
    
    numericInput("K.n", "Carrying capacity", 
                 value = 100, min = 1, step = 1),
    
    numericInput("N0", "Initial number of local populations (N0)", 
                 value = 10, min = 1, step = 1),

    numericInput("r.n", "Growth rate (r)", 
                 value = 1, min = 0, step = .01),
    
    numericInput("alpha", "Per capita assisination by pegamunk (alpha)", 
                 value = 1, min = 0, step = .01),
    
    h3("Pegamunk parameters:"),
    
    numericInput("K.m", "Carrying capacity", 
                 value = 100, min = 1, step = 1),
    
    numericInput("M0", "Initial number of local populations (N0)", 
                 value = 10, min = 1, step = 1),
    
    numericInput("r.m", "Growth rate (r)", 
                 value = 1, min = 0, step = .01),
    
    numericInput("beta", "Per capita murder by squirlicorn (beta)", 
                 value = 1, min = 0, step = .01),
    
    h3("Prepare the stage:"),
    numericInput(inputId = "tmax", 
                 label = "Time:", 
                 value = 2, min = 2, max = 400),
    
    actionButton("go", "May the battle begin!")),
  
  fluidRow( 
    column(width = 7, style='padding:10px', plotOutput("plots", height = '800px', width = '800px')))
)

server <- function(input, output) {
  sim <- eventReactive(input$go,
                       runCompetition(r.n = input$r.n, 
                                      r.m = input$r.m, 
                                      K.n = input$K.n, 
                                      K.m = input$K.m, 
                                      N0 = input$N0, 
                                      M0 = input$M0, 
                                      alpha= input$alpha, 
                                      beta = input$beta, 
                                      tmax = 200))
  output$plots <- renderPlot({
    mysim <- sim()
    plotMN(input$tmax, mysim)
  })
}

shinyApp(ui=ui, server=server)