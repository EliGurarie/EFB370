

plotIso <- function(K1, alpha, K2, beta){
    pars(); par(xaxs="i", pty = "m")
    par(mgp = c(1.75,.75,0))
    plot(0,0, type = "n", xlim = c(0,max(K1,K2/beta)), ylim = c(0,max(K2,K1/alpha)), 
         xlab = expression(N[1]), ylab = expression(N[2]), 
         xaxt = "n", yaxt = "n")
    
    abline(K1/alpha, -1/alpha, add = TRUE, col = "darkorange", lwd = 3)
    points(c(0,K1), c(K1/alpha, 0), pch = 19 , cex = 2, col = "darkorange", xpd = NA)
    
    abline(K2, -beta, add = TRUE, col = "darkblue", lwd = 3)
    points(c(K2/beta, 0), c(0,K2), pch = 19 , cex = 2, col = "darkblue", xpd = NA)
    
    N1star <- -(K1 - alpha*K2)/(beta * alpha - 1)
    N2star <- K2 - beta * N1star
    points(N1star, N2star, pch = 19, col = "red", cex = 2)
    
    axis(side = 1, at = K1, expression(K[1]),
         font = 3, col.axis = "darkorange")
    axis(side = 2, at = K1/alpha, expression(K[1]/alpha),
         font = 3, col.axis = "darkorange", las = 1)
    axis(side = 1, at = K2/beta, expression(K[2]/beta),
         font = 3, col.axis = "darkblue")
    axis(side = 2, at = K2, expression(K[2]),
         font = 3, col.axis = "darkblue", las = 1)
}




plotPredPreyPhase <- function(K1, alpha, K2, beta){
    pars(); par(xaxs="i", pty = "m")
    par(mgp = c(1.75,.75,0))
    plot(0,0, type = "n", xlim = c(0,K2)*1.5, ylim = c(0,K1)*1.5, 
         ylab = "predator (P)", xlab = "prey (V)", 
         xaxt = "n", yaxt = "n")
    
    abline(K1/alpha, -1/alpha, add = TRUE, col = "darkorange", lwd = 3)
    points(c(0,K1), c(K1/alpha, 0), pch = 19 , cex = 2, col = "darkorange", xpd = NA)
    
    abline(K2, -beta, add = TRUE, col = "darkblue", lwd = 3)
    points(c(K2/beta, 0), c(0,K2), pch = 19 , cex = 2, col = "darkblue", xpd = NA)
    
    axis(side = 1, at = K1, expression(gamma/r[p]),
         font = 3, col.axis = "darkorange")
    axis(side = 2, at = K2, expression(sigma/r[v]),
         font = 3, col.axis = "darkblue", las = 1)
    
    
    N1star <- -(K1 - alpha*K2)/(beta * alpha - 1)
    N2star <- K2 - beta * N1star
    points(N1star, N2star, pch = 19, col = "red", cex = 2)
}
