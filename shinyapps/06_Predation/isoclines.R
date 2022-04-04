# predator prey model
rm(list=ls())
library(deSolve)
require(magrittr)
require(shiny)
require(shinyjs)
library(shinyWidgets)

Isoclines <- function(competition, predation, rv, rp, K, xlim, n.arrows){
  
  alpha <- beta <- competition
  gamma <- sigma <- predation
  
  pars <- list(alpha= alpha, beta = beta, 
               gamma = gamma, sigma = sigma, 
               K = K, rv = rv, rp = rp)
  
  getdPdN <- function(v, p, pars){
    dP <- with(pars, rp * p * (1 - p/K - beta * v/K) + sigma * v * p)
    dV <- with(pars, rv * v * (1 - v/K - alpha * p/K) - gamma * v * p)
    dV + 1i*dP
  }

  par(mar = c(3,3,2,2), bty = "l", mgp = c(1.5,.25,0),
      tck = 0.01, cex.lab = 1, cex.axis = 0.8, xaxs="i", yaxs = "i", las = 1)
  plot(0,0, type = "n", col = "purple", 
       ylim = c(0, xlim), xlim = c(0, xlim),
       xlab = "prey", ylab = "predator", asp = 1)

  abline(v = K, col = "darkgrey", lwd = 2, lty = 3)
  abline(h = K, col = "darkgrey", lwd = 2, lty = 3)
  
  if(sigma > 0) abline(h = rv/sigma, lwd = 2, col = "darkblue")
  cols <- c("darkorange", "darkblue")
  
  if(competition == 0 & predation == 0){
    abline(v = K, col = cols[1], lwd = 3)
    abline(h = K, col = cols[2], lwd = 3)
  } else {
    a.prey <- rv*K/(alpha*rv + K*gamma)
    b.prey <- -rv/(alpha*rv + K*gamma)
    a.pred <- K
    b.pred <- sigma*K/rp - beta
    #a.prey <- K/(beta + sigma*K/rp)
    #b.prey <- -1/(beta + sigma*K/rv)
    #a.pred <- K
    #b.pred <- -alpha + gamma*K/rp
    abline(a.prey, b.prey, col =  cols[2], lwd = 3)
    abline(a.pred, b.pred, col =  cols[1], lwd = 3)
  }
  
  z.grid <- outer(xlim*seq(0,1.4,length = n.arrows), 
                  xlim*seq(0,1.4,length = n.arrows), 
                  function(x,y) x + 1i*y) 
  
  arrow <- getdPdN(Re(z.grid), Im(z.grid), pars = pars)/5
  arrow.scaled.x <- (Re(arrow)-min(Re(arrow)))/diff(range(Re(arrow)))
  arrow.scaled.y <- (Im(arrow)-min(Im(arrow)))/diff(range(Im(arrow)))
  
  cols <- rgb((arrow.scaled.x * as.numeric(Re(arrow) > 0))/2 + 
                (arrow.scaled.y * as.numeric(Im(arrow) > 0))/2, 
              (arrow.scaled.x * as.numeric(Re(arrow) < 0))/2 + 
                (arrow.scaled.y * as.numeric(Im(arrow) > 0))/2,  
              (arrow.scaled.y * as.numeric(Im(arrow) > 0))/2 + 
                (arrow.scaled.y * as.numeric(Im(arrow) < 0))/2)

  z.end <- z.grid + arrow
  arrows(Re(z.grid), Im(z.grid), Re(z.end), Im(z.end), 
         length = .05, col =cols, lwd = 2)
}

#Isoclines(competition = .01, 
#          predation = .01, 
#          rp = .1, rv = .1, 
#          K = 1)


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
  tags$style(
    "*, div {
      font-family: Verdana;
    }"
  ),
  setBackgroundColor("Black"),
  h1("The Original Isocline Mega-Mapper", style = "color: #00FF00"),
  sidebarLayout(
    sidebarPanel(
      h2("Sprinkling Colorful Arrows Around Phase Spaces Since 2022!", style = "color: #006432"),
      h3("Processes:"),
      sliderInput("competition", withMathJax("Competition intensity  ($\\alpha = \\beta$)"), 
                   value = 0, min = 0, step = .01, max = 5),
      sliderInput("predation", withMathJax("Predation intensity  ($\\sigma = \\gamma$)"), 
                  value = 0, min = 0, step = .001, max = 2),
      
      fluidRow(
        column(6, 
               sliderInput("rv", withMathJax("Prey growth ($r_v$)"), 
                           value = .1, min = 0, step = .01, max = 2)),
        column(6, 
             sliderInput("rp", withMathJax("Predator growth ($r_p = -q$)"), 
                         value = .1, min = -1, step = .01, max = 1))),
      numericInput("K", withMathJax("Carrying capacity"), value = 1, min = 0, step = 1),
      h3("Visuals:"),
      fluidRow(
        column(6, 
               numericInput("xlim", withMathJax("Plot range"), 
                            value = 1.3, min = 0, step = .01)),
        column(6, 
               numericInput("n.arrows", withMathJax("Number of arrows"), 
                            value = 20, min = 4, step = 1))
     )),
    mainPanel(imageOutput("plots", height = '92%', width = '92%'))
  )
)

server <- function(input, output) {
  output$plots <- renderPlot(
    Isoclines(input$competition, input$predation, input$rv,  input$rp, input$K, input$xlim, input$n.arrows),
    res = 144, width = 792, height = 792)
}

shinyApp(ui=ui, server=server)