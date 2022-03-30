# predator prey model
rm(list=ls())
library(deSolve)
require(magrittr)
require(shiny)
require(shinyjs)
library(shinyWidgets)

runLV <- function(tmax, dt=.01,  
                  parms = c(gamma, sigma, r, q, p0, v0)){
    LV <- function(t, y, parameters) {
        with(as.list(c(y, parameters)), {
            
            dV <- r * V - sigma * P * V
            dP <- gamma * P * V - q * P
            
            return(list(c(dP, dV)))
        })
    }
    yini <- c(P = as.numeric(parms["p0"]), 
              V = as.numeric(parms["v0"]))
    times <- seq(from = 0, to = tmax, by = dt)
    list(out = ode(y = yini, times = times, func = LV, parms = parms),
         parms = parms)
}

# parms = c(gamma=1, sigma=1, r=1, q=1, p0=.5, v0=2)
# a <- runLV(40,.01,parms)
# plotLV(a, 5)

plotLV <- function(sim, t.show){
    par(mar = c(3,3,2,2), bty = "l", mgp = c(1.5,.25,0),
        tck = 0.01, cex.lab = 1.2)
    layout(1:2 %>% t, widths = 2:1)
    cols <- c("darkorange", "darkblue")
    
    out <- sim$out[sim$out[,"time"] < t.show,]
    parms <- sim$parms
    
    plot(out[,"time"],out[,"P"], type = "l", xlab = "time", 
         ylim = c(0, max(out[,2:3])), 
         ylab = "abundance", col = cols[1], lwd = 2)
    lines(out[,"time"],out[,"V"], col = cols[2], lwd = 2)
    legend("topleft", legend = c("predator", "prey"), lwd = 2, col = cols)
    
    plot(out[,"V"],out[,"P"], type = "l", col = "purple", 
         ylim = c(0, max(out[,"P"])),
         xlim = c(0, max(out[,"V"])),
         xlab = "prey", ylab = "predator", lwd = 2)
    points(out[,"V"] %>% tail(1), out[,"P"] %>% tail(1), 
           col = "purple", bg = "pink", pch = 21, cex = 3)
    abline(h = parms["r"]/parms["sigma"], col = cols[2], lty = 3, lwd = 3)
    abline(v = parms["q"]/parms["gamma"], col = cols[1], lty = 3, lwd = 3)
    points(parms["q"]/parms["gamma"], parms["r"]/parms["sigma"],  
           pch = 19, col = "darkred", cex = 2)
}


ui <- fluidPage(
    tags$head(
        tags$link(rel="stylesheet", 
                  href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css", 
                  integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
                  crossorigin="anonymous"),
        HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
        HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
        HTML('
      <script>
        document.addEventListener("DOMContentLoaded", function(){
          renderMathInElement(document.body, {
            delimiters: [{left: "$", right: "$", display: false}]
          });
        })
      </script>')
    ),
    setBackgroundColor("Black"),
    h1("Lotka-Volterra Predator-Prey", style = "color: Yellow"),
    
    sidebarPanel(
        h3("Prey process:"),
        numericInput("v0", "Initial number of prey ($V(0)$)", 
                     value = 5, min = 1, step = 1),
        numericInput("r", "Intrinsic growth rate ($r$)", 
                     value = 1, min = 0, step = .01),
        numericInput("gamma", withMathJax("Predator consumption ($\\gamma$)"), 
                     value = 1, min = 0, step = .01),
        
        h3("Predator process:"),
        numericInput("p0", "Initial number of predators ($P(0)$)", 
                     value = 3, min = 1, step = 1),
        numericInput("q", "Mortality of predators ($N_0$)", 
                     value = 10, min = 1, step = 1),
        numericInput("sigma", "Predator growth from prey consumption ($\\sigma$)", 
                     value = 1, min = 0, step = .01),

        h3("Prepare the stage:"),
        numericInput(inputId = "t", 
                     label = "Time:", 
                     value = .1, min = .1, max = 40, step = .1),
        actionButton("go", "Run process!")),
    
    fluidRow(column(width = 7, style='padding:10px', 
                    plotOutput("plots", height = '400px', width = '1000px')))
)

server <- function(input, output) {
    sim <- eventReactive(input$go, {
        runLV(tmax = 40, dt = 0.01, 
              parms = c(gamma = input$gamma, sigma = input$sigma, 
                        r = input$r, q = input$q, 
                        p0 = input$p0, v0 = input$v0))})
    
    output$plots <- renderPlot({
        mysim <- sim()
        plot(mysim$out[,1:2])
        plotLV(mysim, t.show = input$t)
    })
}

shinyApp(ui=ui, server=server)