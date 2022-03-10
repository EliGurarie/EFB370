# stages <- 1:3
# fecundities <- c(0,1.5,.5)
# mortalities <- c(.5,0,1)
# start <- c(20,0,0)

require(expm)
require(gplots)
require(kableExtra)
require(plyr)
require(magrittr)
require(shinyMatrix)
require(utils)


runPop <- function(lifehistory, tmax){
    
    n.stages <- nrow(lifehistory)
    births <- lifehistory[,"births"]
    survival <- lifehistory[,"survival"]
    
    M <- matrix(0, nrow = n.stages, ncol = n.stages)
    M[1,] <- births
    for(i in 2:n.stages) M[i,i-1] <- survival[i-1]
    row.names(M) <- row.names(lifehistory)
    
    time <- 0:tmax
    pop <- sapply(time, function(t) {(M%^%t) %*% lifehistory[,"start"]})
    row.names(pop) <- row.names(lifehistory) 
    
    return(list(pop = pop, leslie = M))
}

m <-cbind(c(20,0,0),c(.5,1,0),c(0,1.5,.5))
colnames(m) <- c("start", "survival", "births")
row.names(m) <- as.roman(1:3) %>% as.character
 
ui <- fluidPage(
    sidebarPanel(
        h4("Life history parameters:"),
        matrixInput(
            "lifehistory",
            value = m,
            rows = list(extend = TRUE, 
                        names = TRUE, 
                        editableNames = TRUE, 
                        delta = 2,
                        delete = TRUE),
            cols = list(names = TRUE, delta = 1),
            class = "numeric"),
            
        numericInput(inputId = "tmax", 
                     label = "Duration:", 
                     value = 20),
        h3(actionButton("button", "Update Life Table")),
        h3(actionButton("go", "Run Simulation"))),
    
    mainPanel(
        plotOutput("structuredgrowth", width = "100%", height = "350px"),
        h3("Leslie Matrix"),
        tableOutput("lesliematrix")
    )
)

server <- function(input, output) {
    sim <- eventReactive(input$go,
                         runPop(input$lifehistory, 
                                input$tmax))
    
    output$structuredgrowth = renderPlot({
        
        mysim <- sim()
        pop <- mysim$pop
        
        palette(rich.colors(nrow(pop)))
        par(mfrow = c(1,2), tck = 0.01, bty = "l", las= 1, 
            cex.lab = 1.5, mgp = c(1.25,.25,0))
        
        layout(t(1:2), widths = c(2,1))
        
        matplot(t(pop), type = "o", pch = 19, lty = 1, xlab = "time", ylab = "N")
        
        legend("top", col = 1:nrow(pop), ncol = ceiling(nrow(pop)*3/3), 
               lty = 1, legend = row.names(pop), title = "age class")
        
        barplot(pop[,ncol(pop)], col = 1:nrow(pop), 
                main = "stable age distribution") 
        
    })
    
   output$lesliematrix = renderTable(sim()$leslie)
}


shinyApp(ui=ui, server=server)
