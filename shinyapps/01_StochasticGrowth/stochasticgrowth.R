
demographicstochasticity <- function(t, N0, b.prob, d.prob, N.sims,
                             whichfirst = "birth"){
  
  Sims <- matrix(NA, nrow = t+1, ncol = N.sims)
  Sims[1,] <- N0
  
  s.prob <- 1-d.prob # survival probability
  
  for(j in 1:ncol(Sims)){
    for(i in 2:(t+1)){
      N_t <- Sims[i-1,j]
      if(whichfirst == "birth"){
        if(N_t < 1e2){
          N_offspring <- rbinom(1, N_t, b.prob)
          N_survive <- rbinom(1, N_t + N_offspring, s.prob)
        } else {
          N_offspring <- round(rnorm(1, mean = b.prob * N_t, 
                               sd = sqrt(N_t * b.prob * (1-b.prob))))
          N_survive <- round(rnorm(1, mean = s.prob * (N_t + N_offspring), 
                             sd = sqrt((N_t + N_offspring) * s.prob * (1-s.prob))))
        }
        Sims[i,j] <- N_survive
      } else {
        N_survive <- rbinom(1, N_t, s.prob)
        N_offspring <- rbinom(1, N_survive, b.prob)
        Sims[i,j] <- N_survive + N_offspring
      }
    }
  }
  extinctions <- apply(Sims, 2, function(x) min(which(x==0)))
  extinctions <- extinctions[extinctions < Inf]
  pop.median <- apply(Sims,1, median)
  pop.low <- apply(Sims,1, quantile, .25)
  pop.high <- apply(Sims,1, quantile, .95)
  
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
                        demographicstochasticity(
                           t = input$t,
                           N.sims = input$N.sims,
                           N0 = input$N0,
                           b.prob = input$b.prob,
                           d.prob = input$d.prob,
                           whichfirst = input$whichfirst))

    output$timeseries <- renderPlot({
        t <- nrow(sim()$simulation)-1
        par(mfrow = c(1,2), tck = 0.01, bty = "l", las= 1, 
            cex.lab = 1.5, mgp = c(1.5,.25,0))
        layout(rbind(1:2,1:2), width = c(2,1))
        matplot(0:t, sim()$simulation, type = "l", lty = 1, col = rgb(0,0,0,.3), 
                xlab = "Time", ylab = "Population")
        points(sim()$extinctions, rep(0, sim()$n.extinct), col = "green", pch = 4, cex = 2, lwd = 2)
        lines(0:t, sim()$pop.median, col = "red", lwd = 5)
        lines(0:t, sim()$pop.low, col = "darkred", lwd = 2)
        lines(0:t, sim()$pop.high, col = "red", lwd = 2)
    })
    
    output$p.extinction <- renderText({
      paste("Probability of extinction:", sim()$n.extinct/ncol(sim()$simulation))
    })
    
    output$survivalreport <- renderText({
      paste("Final population - mean:", round(mean(sim()$final.pop),3),
            "\r\nFinal population - standard deviation:", round(sd(sim()$final.pop), 3))
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
        numericInput("b.prob", label = withMathJax(helpText("$p_b$- probability of birth")), value = .5, step = .01, min = 0, max = 1),
        numericInput("d.prob", label = withMathJax("$p_d$ - probability of death"), value = 0.5, step = .01, min = 0, max = 1),
        radioButtons("whichfirst", "Which is first? ", c("Birth" = "birth", "Death" = "death")),

        h4("Simulation properties"),
        numericInput("N0", label = "$N_0$ - initial population", value = 10),
        numericInput("t", label = withMathJax("$t_{max}$ - number of years"), value = 20),
        numericInput("N.sims", label = withMathJax("Number of simulations"), value = 100),
        
        h2(actionButton("go", "Run simulation"))),
  
    mainPanel(
        plotOutput("timeseries", height = "400px", width = "1200px"),
        
        h4("Probability of Extinction:"),
        verbatimTextOutput("p.extinction"),
        tags$head(tags$style("#report{color:darkblue; font-size:14px; background: ghostwhite;}")),
        
        h4("Survival Report:"),
        verbatimTextOutput("survivalreport"),
        tags$head(tags$style("#report{color:darkblue; font-size:14px; background: ghostwhite;}")),
        
        h4("Extinction Report:"),
        verbatimTextOutput("extinctionreport"),
        tags$head(tags$style("#report{color:darkblue; font-size:14px; background: ghostwhite;}"))
  )
    
)

# shinyApp(ui=ui, server=server)
