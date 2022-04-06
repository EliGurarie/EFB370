require(deSolve)

runLV <- function(tmax, dt=.01,  
                  parms = c(r, q, gamma, sigma, delta, p0, v0, type)){

  FR1 <- function(v, sigma, delta) {sigma*v} 
  FR2 <- function(v, sigma, delta) {sigma*v / (v + delta)} 
  FR3 <- function(v, sigma, delta) {sigma*v^2 / (v^2 + delta^2)}
  FR4 <- function(v, sigma, delta) {sigma*v^2 / (v^2 + delta^2)}
  
  FunctionalResponse <- get(paste0("FR",type))
  
   PPwithFR <- function(t, y, parameters) {
  
    with(as.list(c(y, parameters)),{
      #if(type == 1) dV <- r * V  - sigma * V * P
      #if(type == 2) dV <- r * V  - sigma * V * P / (V + delta)
      #if(type == 3) dV <- r * V  - sigma * V^2 * P / (V^2 + delta^2)
      dV <- r * V  - FunctionalResponse(V, sigma, delta) * P
      dP <- gamma * P * V - q * P
      return(list(c(dP, dV)))
    })
  } 

  yini <- c(P = as.numeric(parms["p0"]), 
            V = as.numeric(parms["v0"]))
  times <- seq(from = 0, to = tmax, by = dt)
  list(out = ode(y = yini, times = times, func = PPwithFR, parms = parms),
       parms = parms)
}

parms = c(gamma=1, sigma=1, delta = 0.5, r=.5, q=1, p0=1.5, v0=.1, type = 3)
a <- runLV(50, .01,  parms)
plotLV(a,50)

curve(FunctionalResponse(x, parms["sigma"], parms["delta"]), add  = TRUE)


plotLV <- function(sim, t.show -= ){
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
    a <- parms["r"] /parms["gamma"]
    b <-  -parms["r"]/parms["Kv"] * parms["gamma"]
    
    #a.prey <- rv*K/(alpha*rv + K*gamma)
    #b.prey <- -rv/(alpha*rv + K*gamma)
    
    abline(a, b, col =  scales::alpha(cols[2], .4), lwd = 3)
    abline(v = x.iso, col = scales::alpha(cols[1], .4), lwd = 3)
    points(x.iso, a + b*x.iso, pch = 19, col = "darkred", cex = 2)
  } else {
    abline(h = y.iso, col = scales::alpha(cols[2], .4), lwd = 3)
    abline(v = x.iso, col = scales::alpha(cols[1], .4), lwd = 3)
    points(x.iso, y.iso, pch = 19, col = "darkred", cex = 2)
  }
}
