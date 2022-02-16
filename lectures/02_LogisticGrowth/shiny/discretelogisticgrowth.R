library("deSolve")

brusselator <- function(t, y, p) {
    with(as.list(c(y, p)), {
        dX <- k1*A   - k2*B*X    + k3*X^2*Y - k4*X
        dY <- k2*B*X - k3*X^2*Y
        list(c(X=dX, Y=dY))
    })
}


logisticgrowth <- function(t, N0, r, K){
    N <- rep(NA, t+1)
    N[1] <- N0
    Time <- 0:t
    for(i in 2:(t+1))
        N[i] <- N[i-1] + N[i-1] * r * (1 - N[i-1]/K)
    data.frame(Time, N)
}

server <- function(input, output) {
    output$logistic <- renderPlot({
        parms <- c(r=input$r, N0=input$N0, K = input$K, 
                   t = input$t, 
                   breaks = input$histbreaks,
                   histstart = input$histstart)
        out <- logisticgrowth(parms["t"],
                              parms["N0"],
                              parms["r"],
                              parms["K"])
        par(mfrow = c(1,2), tck = 0.01, bty = "l", las= 1, 
            cex.lab = 1.5, mgp = c(1.25,.25,0))
        layout(rbind(1:2,1:2), width = c(2,1))
        plot(out, type = "o", cex = 1.5, pch = 19, col="grey")
        

        hist(out$N[parms["histstart"]:(parms["t"]+1)], 
             breaks = parms["breaks"], 
             col = "grey", xlab = "population size",
             main = "Histogram of values")

    })
}

ui <- fluidPage(
    
    
    sidebarPanel(
        h4("Simulation properties"),
        numericInput("N0", label = "N0 - initial population", value = 10),
        numericInput("t", label = "t - number of years", value = 2),
        
        h4("Parameters"),
        numericInput("r", label = "r - growth rate", value = .5, step = .01),
        numericInput("K", label = "K - carrying capacity", value = 100))
        
        #h4("Plotting"),
       ,
    mainPanel(
        #h3("Logistic growth"),
        plotOutput("logistic", height = "350px"),
        numericInput("histbreaks", label = "histogram breaks", value = 30, step = 1),
        numericInput("histstart", label = "histogram start", value = 1, step = 1)
    )
    
)

#shinyApp(ui=ui, server=server)
