# predator prey model
library(deSolve)
require(magrittr)
require(shiny)
require(shinyjs)
library(shinyWidgets)

runWCM <- function(tmax, dt=.01, parms){
  # parms contains: 
  # c(sigma_m, sigma_c, delta, 
  #   alpha_cm, alpha_mc, 
  #   gamma_c, gamma_m, 
  #   rc, rm, Km, Kc,  w0, m0, c0)
  
  LV <- function(t, y, parameters) {
    with(as.list(c(y, parameters)), {
      
      dW <- gamma_m * M * W + gamma_c * C * W - delta * W
      dM <- rm*(1 - M/Km) * M - sigma_m * W * M
      dC <- rc*(1 - C/Kc) * C - sigma_c * W * C
      
      return(list(c(dW, dM, dC)))
    })
  }
  
  yini <- c(W = as.numeric(parms["w0"]), 
            M = as.numeric(parms["m0"]),
            C = as.numeric(parms["c0"]))
  times <- seq(from = 0, to = tmax, by = dt)
  list(out = ode(y = yini, times = times, func = LV, parms = parms),
       parms = parms)
}



plotWCM <- function(sim, t.show = max(sim$out[,1])){
  par(mar = c(3,3,2,2), bty = "l", mgp = c(1.5,.25,0),
      tck = 0.01, cex.lab = 1.2)
  
  layout(rbind(c(1,1,1), 2:4), heights = c(1.5,1))

  out <- sim$out[sim$out[,"time"] < t.show,]
  
  time <- out[,"time"]
  w <- out[,"W"]
  c <- out[,"C"]
  m <- out[,"M"]
  
  cols <- c("grey", "red", "green")
  outercols <- c("darkgrey", "darkred", "darkgreen")
  
  plot(time, w, type = "l", xlab = "time", 
       ylim = c(0, max(w, c, m, na.rm = TRUE)), 
       ylab = "abundance", col = cols[1], lwd = 2)
  lines(time, m, col = cols[2], lwd = 2)
  lines(time, c, col = cols[3], lwd = 2)
  points(rep(time[length(time)], 3), c(w[length(w)], m[length(m)], c[length(c)]), 
         pch = c(24, 21, 22), cex = 2, col = outercols, bg = cols)
  
  legend("top", legend = c("wolf", "moose", "caribou"), lwd = 2, col = cols, ncol = 3)
  
  plot(c, w, type = "l", col = "purple", 
       ylim = c(0, max(w, na.rm = TRUE)),
       xlim = c(0, max(c, na.rm = TRUE)),
       xlab = "caribou", ylab = "wolf", lwd = 2)
  points(c %>% tail(1), w %>% tail(1), 
         col = "purple", bg = "antiquewhite", pch = 21, cex = 2)
  
  plot(m, w, type = "l", col = "purple", 
       ylim = c(0, max(w, na.rm = TRUE)),
       xlim = c(0, max(m, na.rm = TRUE)),
       xlab = "moose", ylab = "wolf", lwd = 2)
  points(m %>% tail(1), w %>% tail(1), 
         col = "purple", bg = "antiquewhite", pch = 21, cex = 2)
  
   plot(c, m, type = "l", col = "purple", 
       ylim = c(0, max(m, na.rm = TRUE)),
       xlim = c(0, max(c, na.rm = TRUE)),
       xlab = "caribou", ylab = "moose", lwd = 2)
  points(c %>% tail(1), m %>% tail(1), 
         col = "purple", bg = "antiquewhite", pch = 21, cex = 2)
}

parms = c(gamma_m = .003, gamma_c = .0001, delta = .15, w0 = 1000, 
          sigma_m = .003, rm = 1, Km = 10000, m0 = 0, 
          sigma_c = .001, rc = 1, Kc = 10000, c0 = 1000)
#a <- runWCM(200,.01,parms)
#plotWCM(a)

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
  
  setBackgroundColor("DarkBrown"),
  h1("Caribou v. Moose v. Wolf: The Three-way Throwdown", style = "color: LimeGreen"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(column(6, 
                      h3("Caribou:"),
                      numericInput("c0", "Initial pop ($C(0)$)", value = 100, min = 0, step = 1),
                      numericInput("rc", "Growth rate ($r_c$)",  value = 1, min = 0, step = .01),
                      numericInput("sigma_c", "Predator mortality ($\\sigma_c$)", value = .01, min = 0, step = .01),
                      numericInput("Kc", "Carrying capacity ($K_c$)", value = 5000, min = 1, step = 1)
      ),
      column(6, 
             h3("Moose:"),
             numericInput("m0", "Initial pop ($M(0)$)", value = 1, min = 0, step = 1),
             numericInput("rm", "Growth rate ($r_m$)",  value = 1, min = 0, step = .01),
             numericInput("sigma_m", "Predator mortality ($\\sigma_m$)", value = .01, min = 0, step = .01),
             numericInput("Km", "Carrying capacity ($K_m$)", value = 5000, min = 1, step = 1)
      )),
      
      h3("Wolf:"),
      numericInput("w0", "Initial pop ($W(0)$)", value = 20, min = 0, step = 1),
      fluidRow(
        column(6, numericInput("gamma_c", "Growth from eating Prey 1 ($\\gamma_c$)", value = .001, min = 0, step = .005)),
        column(6, numericInput("gamma_m", "Growth from eating Prey 2 ($\\gamma_m$)", value = .003, min = 0, step = .005))),
      numericInput("delta", "Death rate", value = .2, min = 0, step = .01),
      
      h3("Prepare the stage:"),
      sliderInput(inputId = "t", label = "Time:", value = 10, min = 1, max = 200, step = 1),
      actionButton("go", "Run process!")),
    
    mainPanel(
      plotOutput("plots", height = '700px', width = '1000px'),
      h3("Final Counts:", style = "color: darkred"),
      textOutput("finalcounts"),
      tags$head(tags$style("#finalcounts{
                                 font-size: 20px;
                                 }"
      ))
    )
))



server <- function(input, output) {
  sim <- eventReactive(input$go, {
    runWCM(tmax = 200, dt = 0.1, 
           parms = c(sigma_m = input$sigma_m, sigma_c = input$sigma_c, delta = input$delta, w0 = input$w0,
                     gamma_m = input$gamma_m, rm = input$rm, Km= input$Km, m0 = input$m0, 
                     gamma_c = input$gamma_c, rc = input$rc, Kc = input$Kc, c0 = input$c0)) 
  }
  )
  
  output$plots <- renderPlot({
    #plot(sim()$out[,1], type = "o", main = sim()$out)
    plotWCM(sim(), t.show = input$t)
  }, res = 144)
  
  output$finalcounts <- renderText({
    myout <- sim()$out
    results <- myout[myout[,"time"] <= input$t,] %>% tail(1)
    paste0("Wolves: ", round(results[,"W"], 1), 
           ";  Caribou: ", round(results[,"C"], 1),
           ";  Moose: ", round(results[,"M"], 1))
  })
}

shinyApp(ui=ui, server=server)