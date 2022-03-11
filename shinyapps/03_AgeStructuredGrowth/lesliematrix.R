require(expm)
require(gplots)
require(plyr)
require(magrittr)
require(shinyMatrix)
require(utils)
require(shinyjs)
library(shinyWidgets)
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
  
  useShinyjs(),
  setBackgroundColor("antiquewhite"),
  tags$style(
    "*, div {
      font-family: Consolas;
    }"
  ),
  
  h1("Matrix Modeller 5000"),
  
  sidebarPanel(
    h3("Enter Leslie matrix:"),
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
    
    h3("Enter initial population:"),
    h5("Must be the same length as the Leslie matrix has columns:"),
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
    actionButton("eigen", "Compute eigenvalues"),
    actionButton("go", "Run simulation")),
    
 h3(actionButton("resetAll", "Clear everything")),
    
  mainPanel(fluidRow( 
      column(width = 12, style='padding:10px', grVizOutput("diagram", height = "200px")), 
      column(width = 12, style='padding:10px', align = "center", plotOutput("eigenvalues", height = "100px", width = "80%")),
      column(width = 12, style='padding:10px', plotOutput("structuredgrowthplot", height = "350px"))
  ))
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
    eigenvector <- Re(eigen(M)$vectors[,1])
    eigenvector <- eigenvector/sum(eigenvector)
    eigenvectors <- paste0("{",paste(round(eigenvector, 3), collapse = ", "),"}")
    
    par(mar = c(0,0,4,0), bty = "n", mfrow = c(1,2), cex.main = 2)
    plot(0, 0, type = "l", xaxt = "n", yaxt = "n",  main = expression("eigenvalue "~lambda~":"))
    text(0, 0, round(Re(eigenvalue),3), cex = 2)
    plot(0, 0, type = "l", xaxt = "n", yaxt = "n",  main = expression("eigenvector "~N^"*"~":"))
    text(0, 0, eigenvectors, cex = 2*sqrt(3/length(eigenvector)))
  })
  
  simplot <-  eventReactive(input$go, {
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

  
  observeEvent(input$go,{
    output$structuredgrowthplot = renderPlot({simplot()})
  })
  
  observeEvent(input$drawdiagram,{
    output$diagram = renderGrViz({diagram()})
  })
  
  observeEvent(input$eigen,{
    output$eigenvalues = renderPlot({eigenplot()})
  })

  observeEvent(input$resetAll, {
    output$diagram <- renderGrViz({})
    output$eigenvalues <- renderPlot({})
    output$structuredgrowthplot <- renderPlot({})
  })
}


# shinyApp(ui=ui, server=server)



