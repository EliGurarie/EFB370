p <- read.csv("Paramecium.csv")
N.logistic <- function(x, N0, K, r0) K/(1 + ((K - N0)/N0)*exp(-r0*x))

plot(Nt~Day, data = p)
# try different values!
N0 <- 10
K <- 150
r0 <- 1.5
plot(Nt~Day, data = p)
curve(N.logistic(x, N0, K, r0), add = TRUE, col = 2)
