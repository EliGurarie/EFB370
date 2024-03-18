


islandmainland <- function(axis = TRUE, ...){
    c <- e <- 1; pstar <- c/(e+c)
    cols <- c("darkgreen","darkred")
    curve(c-c*x, xlim = c(0,1), ylim = c(0,c),
          lwd = 2, col = cols[1], 
          yaxt = "n", ...)
    abline(0, e, lwd = 2, col = cols[2])
    abline(v = pstar, col = "darkgrey", lwd = 2, lty = 3)
    if(axis){
        axis(1, at = c(0, pstar, 1), 
             labels = c(0, expression(f^'*'), "1"))
        text(c(.15, .85), c(.6,.6), 
             c(expression(c == p[c](1-f)), expression(e == p[e]~f)), 
             col = cols)}
}

internalcolonization <- function(axis = TRUE, ...){
    c <- .4; e <- .1; pstar <- 1-e/c
    cols <- c("darkgreen","darkred")
    curve(c*x*(1 - x), xlim = c(0,1), ylim = c*c(0,.3),
          main = "", lwd = 2, col = cols[1], 
          yaxt = "n", xaxt = "n", ...)
    curve(e*x*(1 - x), add = TRUE, lwd = 2, col = cols[1])
    abline(0, e, lwd = 2, col = cols[2])
    abline(v = c(0,pstar), col = "darkgrey", lwd = 2, lty = 3)
    if(axis){
        axis(1, at = c(0, pstar, 1), 
             labels = c(0, expression(f^'*'), "1"))
        text(c(.5,.5, pstar), c(e/4, .9*c/4, e*pstar), 
             c(expression(low~p[c]),
               expression(c == p[c]*f(1-f)),
               expression(e == p[e]*f)), 
             pos = c(1,1,4),  col = cols[c(1,1,2)])
    }
}

rescueeffect <- function(label = TRUE, ...){
    c <- .2; e <- .6
    e2 <- .2
    cols <- c("darkgreen","darkred")
    pstar = c/e
    curve(c - c*x, xlim = c(0,1), ylim = c*c(0,1),
          lwd = 2, col = cols[1],
          yaxt = "n", ...)
    curve(e*x*(1 - x), add = TRUE, lwd = 2, col = cols[2])
    curve(e2*x*(1 - x), add = TRUE, lwd = 2, col = cols[2])
    if(label)
    text(c(.7,.5), c(c/2, e/12), 
         c(expression(c == p[c](1-f)),
           expression(e == p[e]*f(1-f))), 
         pos = c(3,3),  col = cols)
    abline(v = c(pstar,1), col = "darkgrey", lwd = 2, lty = 3)
}

both <- function(label = TRUE, ...){
    c <- .3; e <- .2
    c2 <- .1
    cols <- c("darkgreen","darkred")
    curve(c*x*(1 - x), xlim = c(0,1), ylim = c*c(0,.3),
          main = "", lwd = 2, col = cols[1], 
          yaxt = "n", ylab = "colonization / extinction rate",
          xlab = "probability of incidence")
    curve(e*x*(1 - x), add = TRUE, lwd = 2, col = cols[2])
    curve(c2*x*(1 - x), add = TRUE, lwd = 2, col = cols[1])
    if(label)
    text(rep(.5,3), c(c, e, c2)/4, 
         c(expression(c == p[c]*f(1-f)),
           expression(e == p[e]*f(1-f)), 
           expression(p_c < p_e)), 
         pos = c(3,3,1),  col = cols)
    abline(v = c(0,1), col = "darkgrey", lwd = 2, lty = 3)
}