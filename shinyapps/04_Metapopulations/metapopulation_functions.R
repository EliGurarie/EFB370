runMetapopulation <- function(K = 20, N0 = 10, 
                              e = .1, c = .2, 
                              tmax = 1000,
                              colonization = "external",
                              extinction = "constant"){
  metapop.matrix <- matrix(ncol = K, nrow = tmax)
  metapop.matrix[1,] <- sample(c(rep(1, N0), rep(0, K-N0))) 
  
  for(i in 2:tmax){
    N.then <-  metapop.matrix[i-1,]
    f <- sum(N.then / K)
    if(extinction == "constant")
      metapop.matrix[i,N.then == 1] <- 1-rbinom(sum(N.then), size = 1, prob = e) else 
        metapop.matrix[i,N.then == 1] <- 1-rbinom(sum(N.then), size = 1, prob = e*(1-f))
    if(colonization == "external")
      metapop.matrix[i,N.then == 0] <- rbinom(sum(N.then == 0), size = 1, prob = c) else 
        metapop.matrix[i,N.then == 0] <- rbinom(sum(N.then == 0), size = 1, prob = c*f)
  }
  
  return(list(mp = metapop.matrix, 
              Z =  runif(K) + 1i*runif(K),
              pars =  c(e = e, c = c),
              model = c(extinction, colonization)))
}


plotMP_theory <- function(sim){
  e.constant <- sim$model[1] == "constant"
  c.external <- sim$model[2] == "external"
  K <- ncol(sim$mp)
  e <- sim$pars["e"]
  c <- sim$pars["c"]
  plot(c(0,K), c(0,max(e,c)), type = "n", xlab = "Number of patches", ylab = "Rate")
  
  if(e.constant & c.external){
    abline(c, -c/K, col = "darkgreen", lwd = 2)
    abline(0, e/K, col = "darkred", lwd = 2)
    abline(v = K * c/(c + e), col = "darkgrey", lty = 3, lwd = 2)
  }
  if(!e.constant & c.external){
    abline(c, -c/K, col = "darkgreen", lwd = 2)
    curve(e*(x/K)*(1-x/K), col = "darkred", lwd = 2, add = TRUE)
    abline(v = K * min(1,c/e), col = "darkgrey", lty = 3, lwd = 2)
  }
  if(e.constant & !c.external){
    curve(c*(x/K)*(1-x/K), col = "darkgreen", lwd = 2, add = TRUE)
    abline(0, e/K, col = "darkred", lwd = 2)
    abline(v = K * max(0,1-e/c), col = "darkgrey", lty = 3, lwd = 2)
  }
  if(!e.constant & !c.external){
    curve(c*(x/K)*(1-x/K), col = "darkgreen", lwd = 2, add = TRUE)
    curve(e*(x/K)*(1-x/K), col = "darkred", lwd = 2, add = TRUE)
    abline(v =c(0,K), col = "darkgrey", lty = 3, lwd = 2)
  } 
}
  
plotMP_timeseries <- function(sim){
  mp <- sim$mp
  e <- sim$pars["e"]
  c <- sim$pars["c"]
  
  rowSums(mp) %>% plot(type = "o", ylim = c(0,ncol(mp)), pch = 21, bg = "grey", 
                       xlab = "Time", ylab = "Occupied sites", col = "darkgrey")
  K <- ncol(mp)
  
  ex <- sim$model[1] == "constant"
  cz <- sim$model[2] == "external"
  
  n.hat <- K*ifelse(ex & cz,  c/(c+e), 
                  ifelse(!ex & cz, min(1,c/e), 
                         ifelse(ex & !cz, max(0, (1-e/c)), 
                                NA)))
  
  linecols <- c("orange", "darkblue","darkred")
  abline(h = c(0, n.hat, ncol(mp)), col = linecols, lwd = 2)
  text(rep((nrow(mp)+1)/2, 3), c(0, n.hat, ncol(mp)), c("min", "theory", "max"), 
       col = linecols, pos = 3, adj = 1, xpd = NA, font = 2)
}
