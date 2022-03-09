
environmentalstochasticity <- function(T.max, N0, r.mean, r.sd, N.sims){
    
    Sims <- matrix(0, nrow = T.max, ncol = N.sims)
    Sims[1,] <- N0
    
    i <- 2
    while(sum(Sims[i-1,]) > 0 & i <= T.max ){
        Sims[i,] <- round(Sims[i-1,] * exp(rnorm(N.sims, r.mean, r.sd)))
        i <- i+1
    }    
    
    extinctions <- apply(Sims, 2, function(x) min(which(x==0)))
    extinctions <- extinctions[extinctions < Inf]
    pop.median <- apply(Sims,1, median)
    pop.low <- apply(Sims,1, quantile, .25)
    pop.high <- apply(Sims,1, quantile, .50)
    
    list(simulation = Sims, 
         n.extinct = length(extinctions), 
         pop.median = pop.median, 
         pop.low = pop.low, 
         pop.high = pop.high, 
         extinctions = extinctions,
         final.pop = Sims[nrow(Sims),])
}

server <- function(input, output) {
    sim <- eventReactive(input$go,
                         environmentalstochasticity(
                             T.max = input$T.max,
                             N.sims = input$N.sims,
                             N0 = input$N0,
                             r.mean = input$r.mean,
                             r.sd = input$r.sd))
    
    output$timeseries <- renderPlot({
        t <- nrow(sim()$simulation)-1
        par(mfrow = c(1,2), tck = 0.01, bty = "l", 
            cex.lab = 1.5, mgp = c(2,.25,0))
        
        matplot(0:t, sim()$simulation, type = "l", lty = 1, col = rgb(0,0,0,.3), 
                xlab = "Time", ylab = "Population")
        points(sim()$extinctions, rep(0, sim()$n.extinct), col = "green", pch = 4, cex = 2, lwd = 2)
        lines(0:t, sim()$pop.median, col = "red", lwd = 5)
        lines(0:t, sim()$pop.low, col = "darkred", lwd = 2)
        lines(0:t, sim()$pop.high, col = "red", lwd = 2)
        
        
        matplot(0:t, sim()$simulation, type = "l", lty = 1, col = rgb(0,0,0,.3), 
                xlab = "Time", ylab = "Population (log-scale)", log = "y")
        points(sim()$extinctions, rep(1, sim()$n.extinct), col = "green", pch = 4, cex = 2, lwd = 2)
        lines(0:t, sim()$pop.median, col = "red", lwd = 5)
        lines(0:t, sim()$pop.low, col = "darkred", lwd = 2)
        lines(0:t, sim()$pop.high, col = "red", lwd = 2)
    })
    
    output$p.extinction <- renderText({
        paste("Probability of extinction:", sim()$n.extinct/ncol(sim()$simulation))
    })
    
    output$survivalreport <- renderText({
        paste("Final population - median:", round(median(sim()$final.pop)),
              "\r\nFinal population - inter-quartile range:", 
              paste(quantile(sim()$final.pop, c(.25,.75)), collapse = "-"))
    })
    
    output$extinctionreport <- renderText({
        paste("Time to extinction - mean:", round(mean(sim()$extinctions),3),
              "\r\nTime to extinction - standard deviation:", round(sd(sim()$extinction), 3))
    })
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
    sidebarPanel(
        h4("Parameters"),
        numericInput("r.mean", label = withMathJax(helpText("$ \\mu_r$ - mean growth rate")), value = .5, step = .01, min = 0, max = 1),
        numericInput("r.sd", label = withMathJax("$\\sigma_r$ - standard deviation of growth rate"), value = 0.5, step = .01, min = 0, max = 1),
        
        h4("Simulation properties"),
        numericInput("N0", label = "$N_0$ - initial population", value = 10),
        numericInput("T.max", label = withMathJax("$t_{max}$ - maximum duration of simulation"), value = 20),
        numericInput("N.sims", label = withMathJax("Number of simulations"), value = 100),
        
        h2(actionButton("go", "Run simulation"))),
    
    mainPanel(
        plotOutput("timeseries", height = "400px", width = "100%"),
        
        h4("Summary Report:"),
        verbatimTextOutput("p.extinction"),
        tags$head(tags$style("#report{color:darkblue; font-size:14px; background: ghostwhite;}")),
        verbatimTextOutput("survivalreport"),
        tags$head(tags$style("#report{color:darkblue; font-size:14px; background: ghostwhite;}")),
        verbatimTextOutput("extinctionreport"),
        tags$head(tags$style("#report{color:darkblue; font-size:14px; background: ghostwhite;}"))
    )
    
)
#shinyApp(ui=ui, server=server)
