p <- read.csv("paramecium")
N.logistic <- function(x, N0, K, r0) K/(1 + ((K - N0)/N0)*exp(-r0*x))

nls(Nt ~ N.logistic(Day, N0, K, r0), data = p, 
    start = list(N0 = 1, K = 200, r0 = 0.75))
