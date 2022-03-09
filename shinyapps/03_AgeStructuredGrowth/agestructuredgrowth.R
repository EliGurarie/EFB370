# stages <- 1:3
# fecundities <- c(0,1.5,.5)
# mortalities <- c(.5,0,1)
# start <- c(20,0,0)

require(expm)
require(gplots)
require(kableExtra)
require(plyr)
require(magrittr)


runPop <- function(n.stages, births, survival, init, tmax){
    
    births <- as.numeric(unlist(strsplit(births,",")))
    survival <- as.numeric(unlist(strsplit(survival,",")))
    init <- as.numeric(unlist(strsplit(init,",")))
    
    M <- matrix(0, nrow = n.stages, ncol = n.stages)
    M[1,] <- births
    for(i in 2:n.stages) M[i,i-1] <- survival[i-1]
        
    time <- 1:tmax
    pop <- sapply(time, function(t) {(M%^%t) %*% init})
    return(pop)
}

getLifeHistoryTable <- function(n.stages, births, survival){
    births <- as.numeric(unlist(strsplit(births,",")))
    survival <- as.numeric(unlist(strsplit(survival,",")))
    data.frame(stage = 1:n.stages, births, survival) %>% mutate(
        cum.survival = cumprod(survival),
        
    )
}

ui <- fluidPage(
    sidebarPanel(
        h4("Life history parameters:"),
        numericInput("n.stages", label = "Number of stages", value = 3, min = 1, step = 1),
        textInput('births', 'Fecundity at stage: (separate by comma)', '0,1.5,.5'),
        textInput('survival', 'Survival at stage: ', '.5,1,0'),
        textInput('init', 'Initial distribution: ', '20,0,0'),
        numericInput(inputId = "tmax", 
                     label = "Duration:", 
                     value = 20),
        h2(actionButton("go", "Run simulation"))),
    
    mainPanel(
        plotOutput("structuredgrowth", width = "100%", height = "350px"),
        h3("Life History Table"),
        tableOutput("lifehistorytable")
    )
)

server <- function(input, output) {
    sim <- eventReactive(input$go,
                         runPop(input$n.stages, 
                                input$births, 
                                input$survival, 
                                input$init, 
                                input$tmax))
    
    output$structuredgrowth = renderPlot({
        
        births <- as.numeric(unlist(strsplit(input$births,",")))
        survival <- as.numeric(unlist(strsplit(input$survival,",")))
        init <- as.numeric(unlist(strsplit(input$init,",")))
        n.stages <- input$n.stages
        
        partest <- c(length(input$births), length(survival), length(init), n.stages)
        validate(need(length(births) == n.stages, 
                      paste("Birth vector must have", n.stages, "elements.")),
                 need(length(survival) >= n.stages-1, 
                      paste("Survival vector must have at least", n.stages-1, "elements.")),
                 need(length(init) == n.stages, 
                      paste("Initial age stage distribution must have", n.stages, "elements.")))
                 
            
        # !all(outer(partest, partest, Vectorize(identical))),
        
        
        # plotting
        pop <- sim()
        
        palette(rich.colors(nrow(pop)))
        par(mfrow = c(1,2), tck = 0.01, bty = "l", las= 1, 
            cex.lab = 1.5, mgp = c(1.25,.25,0))
        
        layout(t(1:2), widths = c(2,1))
        
        matplot(t(pop), type = "o", pch = 19, lty = 1, xlab = "time", ylab = "N")
        legend("top", col = 1:n.stages, ncol = ceiling(n.stages*2/3), 
               lty = 1, legend = 1:n.stages, title = "age class")
        
        barplot(pop[,ncol(pop)], col = 1:n.stages, 
                main = "stable age distribution") 
        
    })
    
    output$lifehistorytable = renderTable(
        getLifeHistoryTable(input$n.stages, 
                            input$births, 
                            input$survival))
}


#shinyApp(ui=ui, server=server)
