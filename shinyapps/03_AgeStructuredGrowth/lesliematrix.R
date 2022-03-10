# stages <- 1:3
# fecundities <- c(0,1.5,.5)
# mortalities <- c(.5,0,1)
# start <- c(20,0,0)

require(expm)
require(gplots)
require(Rage)
require(plyr)
require(magrittr)
require(shinyMatrix)
require(utils)
library(DiagrammeR)
source("leslie_diagram.R")

runPop <- function(M, N0, tmax){
  
  time <- 0:tmax
  pop <- sapply(time, function(t) {(M%^%t) %*% N0})
  row.names(pop) <- row.names(M) 
  
  return(pop)
}

M <- rbind(c(0,1.5,.5),c(.5,0,0),c(0,1,0))
N0 <- c(20,0,0) %>% t
  
colnames(M) <- colnames(N0) <- row.names(M) <- as.roman(1:3) %>% as.character

ui <- fluidPage(
  sidebarPanel(
    h4("Leslie Matrix:"),
    h5("Must be square, i.e. equal number of rows and columns:"),
    matrixInput(
      "lesliematrix",
      value = M,
      inputClass = "numeric",
      rows = list(extend = TRUE, 
                  delete = TRUE,
                  copy = TRUE,
                  names = TRUE,
                  editableNames = TRUE),
      cols = list(extend = TRUE, 
                  paste = TRUE, 
                  delete = TRUE, 
                  names = TRUE),
      class = "numeric"),
    
    h4("Initial population:"),
    h5("Must have same number of rows as Leslie matrix"),
    matrixInput(
      "N0",
      value = N0,
      inputClass = "numeric",
      rows = list(extend = FALSE, n = 1),
      cols = list(extend = TRUE, names = TRUE, delete = TRUE),
      class = "numeric"),
    
    numericInput(inputId = "tmax", 
                 label = "Duration:", 
                 value = 20),
    # h3(actionButton("button", "Update Life Table")),
    actionButton("drawdiagram", "Draw diagram"),
    actionButton("go", "Run simulation"),
    actionButton("eigen", "Compute eigenvalues"),),
    
  mainPanel(
    grVizOutput("diagram", height = "250px", width = "50%"),
    plotOutput("eigenvalues", width = "50%", height = "100px"),
    plotOutput("structuredgrowthplot", width = "100%", height = "350px"),
  )
)

server <- function(input, output) {
  
  sim <- eventReactive(input$go,
                       runPop(input$lesliematrix,
                              t(input$N0),
                              input$tmax))
  
  diagram <- eventReactive(input$drawdiagram,
    leslie_diagram(input$lesliematrix, stages = row.names(input$lesliematrix),
                   height = 300)
  )

  eigenplot <- eventReactive(input$eigen,{
    M <- input$lesliematrix
    eigenvalue <- eigen(M)$values[1]
    eigenvector <- eigen(M)$vectors[,1]
    eigenvector <- eigenvector/sum(eigenvector)
    eigenvectors <- paste0("{",paste(round(eigenvector, 3), collapse = ", "),"}")
    
    par(mar = c(0,0,2,0), bty = "n", mfrow = c(1,2), cex.main = 2)
    plot(0, 0, type = "l", xaxt = "n", yaxt = "n",  main = expression("eigenvalue "~lambda~":"))
    text(0, 0, round(eigenvalue,3), cex = 2)
    plot(0, 0, type = "l", xaxt = "n", yaxt = "n",  main = expression("eigenvector "~N^"*"~":"))
    text(0, 0, eigenvectors, cex = 2)
  })
  
  output$structuredgrowthplot = renderPlot({
    
    pop <- sim()
    
    palette(rich.colors(nrow(pop)))
    par(mfrow = c(1,2), tck = 0.01, bty = "l", las= 1, 
        cex.lab = 1.5, mgp = c(1.25,.25,0))
    
    layout(t(1:2), widths = c(2,1))
    
    matplot(t(pop), type = "o", pch = 19, lty = 1, xlab = "time", ylab = "N")
    
    legend("top", col = 1:nrow(pop), ncol = ceiling(nrow(pop)*3/3), 
           lty = 1, legend = row.names(pop), title = "age class")
    
    barplot(pop[,ncol(pop)], col = 1:nrow(pop), 
            main = "final distribution") 
  })
  
  output$diagram = renderGrViz({
    diagram()
  })
  
  output$eigenvalues = renderPlot({
    eigenplot()
  })
}


shinyApp(ui=ui, server=server)



