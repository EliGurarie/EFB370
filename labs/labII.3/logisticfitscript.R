# first load the paramecium.csv data
N.logistic <- function(x, N0, K, r0) K/(1 + ((K - N0)/N0)*exp(-r0*x))

# try different values!
N0 <- 10
K <- 150
r0 <- 1

plot(Nt ~ Day, data = paramecium)
curve(N.logistic(x, N0, K, r0), add = TRUE, col = 2, lwd = 2)