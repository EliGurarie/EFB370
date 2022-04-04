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
    LVK <- function(t, y, parameters) {
      with(as.list(c(y, parameters)), {
        
        dV <- r * V*(1 - V/Kv) - sigma * P * V
        dP <- gamma * P * V - q * P
        
        return(list(c(dP, dV)))
      })
    }
    
    if(is.na(parms["Kv"])) FUN <- LV else FUN <- LVK
    
    yini <- c(P = as.numeric(parms["p0"]), 
              V = as.numeric(parms["v0"]))
    times <- seq(from = 0, to = tmax, by = dt)
    list(out = ode(y = yini, times = times, func = FUN, parms = parms),
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
    
    x.iso <- parms["q"]/parms["gamma"]
    y.iso <-  parms["r"]/parms["sigma"]
    
    if(!is.na(parms["Kv"])){
      a <- 1/parms["gamma"]
      b <-  -(1/(parms["Kv"] * parms["gamma"]))
      abline(a, b, col =  scales::alpha(cols[2], .4), lwd = 3)
      abline(v = x.iso, col = scales::alpha(cols[1], .4), lwd = 3)
      points(x.iso, a + b*x.iso, pch = 19, col = "darkred", cex = 2)
    } else {
      abline(h = y.iso, col = scales::alpha(cols[2], .4), lwd = 3)
      abline(v = x.iso, col = scales::alpha(cols[1], .4), lwd = 3)
      points(x.iso, y.iso, pch = 19, col = "darkred", cex = 2)
  }
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
    sidebarLayout(
      sidebarPanel(
        fluidRow(column(6, 
                 h3("Prey process:"),
                 numericInput("v0", "Initial number  ($V(0)$)", value = 5, min = 1, step = 1),
                 numericInput("r", "Growth rate ($r$)",  value = 1, min = 0, step = .01),
                 numericInput("gamma", withMathJax("Mortality from predator consumption ($\\gamma$)"), value = 1, min = 0, step = .01),
                 numericInput("Kv", withMathJax("Carrying capacity ($K_v$)"), value = NA, min = 1, step = 1)
                 ),
                 column(6, 
                 h3("Predator process:"),
                 numericInput("p0", "Initial number ($P(0)$)", value = 3, min = 1, step = 1),
                 numericInput("q", "Mortality ($q$)", value = 10, min = 1, step = 1),
                 numericInput("sigma", "Predator growth from prey consumption ($\\sigma$)", value = 1, min = 0, step = .01)
          )),
      h3("Prepare the stage:"),
      sliderInput(inputId = "t", label = "Time:", value = 10, min = .1, max = 40, step = .1),
      actionButton("go", "Run process!")),
     
      mainPanel(plotOutput("plots", height = '400px', width = '1000px'))
    )
)

server <- function(input, output) {
    sim <- eventReactive(input$go, {
        runLV(tmax = 40, dt = 0.01, 
              parms = c(gamma = input$gamma, sigma = input$sigma, 
                        r = input$r, q = input$q, 
                        Kv = input$Kv,
                        p0 = input$p0, v0 = input$v0))})
    
    output$plots <- renderPlot({
        mysim <- sim()
        plot(mysim$out[,1:2])
        plotLV(mysim, t.show = input$t)
    })
}

shinyApp(ui=ui, server=server)