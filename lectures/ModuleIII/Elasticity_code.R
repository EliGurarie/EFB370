f <- c(L = 0, S = 1.5, E = 0.5)


M <- rbind(c( 0, 1.5, .5),
       c(.5,   0,  0),
       c( 0,   1,  0))
row.names(M) <- colnames(M) <- c("L","S","E")

lambda <- eigen(M)$values[1]

# change sophomore fecundity

f_s_test <- seq(1,2,.1)
lambdas <- rep(NA, length(f_s_test))

for(i in 1:length(f_s_test)){
  M_modified <- M
  M_modified["L","S"] <- f_s_test[i]
  lambdas[i] <- eigen(M_modified)$values[1]
}

plot(f_s_test, lambdas, type = "o", pch = 19)
f_s <- M["L","S"]

E <- M*0


for(i in 1:3) for(j in 1:3){
  m_ij <- M[i,j]
  if(m_ij > 0){
    m_test <- seq(m_ij-.2, m_ij+.2,length = 10)
    lambdas <- rep(NA, length(m_test))
    for(m in 1:length(m_test)){
      M_modified <- M
      M_modified[i,j] <- m_test[m]
      lambdas[m] <- Re(eigen(M_modified)$values[1])
    }
    E[i,j] <- m_ij/lambda * lm(Re(lambdas) ~ m_test)$coef[2]
  } 
}

E

