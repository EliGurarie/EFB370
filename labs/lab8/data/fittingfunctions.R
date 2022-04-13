fitLogistic <- function(data, y = "N", time = "Day", N0 = 1, K = 200, r0 = 0.75){
    
    N.logistic <- function(x, N0, K, r0) K/(1 + ((K - N0)/N0)*exp(-r0*x))
    
    Y <- with(data, get(y))
    X <- with(data, get(time))
    
    myfit <- nls(Y ~ N.logistic(X, N0, K, r0),  
        start = list(N0 = N0, K = K, r0 = r0))
    
    summary(myfit)
}

linesLogistic <- function(fit, ...){
    N.logistic <- function(x, N0, K, r0) K/(1 + ((K - N0)/N0)*exp(-r0*x))
    curve(N.logistic(x, N0 =  fit$coefficients[1,1],
                     K= fit$coefficients[2,1],
                     r0= fit$coefficients[3,1]), add = TRUE, ...)
}


simulateCompetition <- function(tmax, tstart = 0, pars, dt = .01){
    Day <- seq(tstart, tmax, dt)
    N2 <- N1 <- rep(Day*0)
    with(pars, {
        N1[1] <- n1
        N2[1] <- n2
        for(i in 2:length(N1)){
            dN1 <- r1 * N1[i-1]*(1 - N1[i-1]/K1 - alpha * N2[i-1]/K1) * dt
            dN2 <- r2 * N2[i-1]*(1 - N2[i-1]/K2 - beta * N1[i-1]/K2) * dt
            N1[i] <- N1[i-1] + dN1
            N2[i] <- N2[i-1] + dN2
        }
        return(data.frame(Day, N1, N2))
    })
}    

simPlot <- function(fit, data = NULL, species1 = "N1", species2 = "N2"){
    cols <- c("darkorange", "darkblue")
    if(!is.null(data))
        matplot(data$Day, data[,c(species1,species2)], type = "o", pch = 19, lty = 1, 
                col = cols,
                xlim = c(0,max(fit$Day))) else
    plot(fit$Day, fit$N1, ylim = range(fit$N1, fit$N2)*1.1, type = "n")
    
    lines(fit$Day, fit$N1, pch = 19, lty = 1, 
          col = cols[1], lwd = 2, type = "l")
    lines(fit$Day, fit$N2, pch = 19, lty = 1, 
          col = cols[2], lwd = 2, type = "l")
    
    legend("topright", legend = c(species1, species2), lwd = 3, col = cols)
}

phasePlot <- function(fit, data, pars, species1, species2){
    cols <- c("darkorange", "darkblue")
    plot(data[,c(species1,species2)], type = "l", col = "grey", 
         xlim = c(0, pars$K1*1.2), ylim = c(0,pars$K2*1.2), asp = 1)
    n <- nrow(data)
    points(data[,c(species1,species2)], col = grey( (n:1)/n*.5 + .3), pch = 19)
    lines(fit$N1, fit$N2, cex = 0.5, lwd = 4, col = "purple")
    
    abline(h = 0, lwd = 3, col = "grey", lty = 3)
    abline(v = 0, lwd = 3, col = "grey", lty = 3)
    abline(pars$K1/pars$alpha, -1/pars$alpha, col = cols[1], lwd = 2)
    abline(pars$K2, -pars$beta, col = cols[2], lwd = 2)
    points(fit$N1[nrow(fit)], 
           fit$N2[nrow(fit)], pch = 4, cex = 3, lwd = 2)
    
    legend("topright", legend = c(species1,species2), lwd = 3, col = cols)
}


getR2 <- function(fit, data, species1, species2){
    
    #var_total <- apply(data[,c(species1,species2)], 2, var)
    SST <- apply(data[,c(species1,species2)], 2, function(x) 
        sum(     (x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
    
    data_sim <- merge(data, fit, by = "Day")[,-1]
    residuals <- cbind(N1 = data_sim[,"N1"] - data_sim[,species1],
                       N2 = data_sim[,"N2"] - data_sim[,species2])
    rSS <- apply(residuals, 2, function(x) sum(x^2))
    
    1-rSS/SST
    
     # var_residual <- apply(residuals, 2, var)
    #    1 - var_residual/var_total
    1-rSS/SST
}

