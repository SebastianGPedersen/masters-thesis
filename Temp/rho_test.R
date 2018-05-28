setwd(Sys.getenv("masters-thesis"))
source('Estimation/rho.R')

simAR<-function(n, rho, sd){
  X<-rep(NA,n)
  X[1]<-rnorm(1, 0, sd/(1-rho))
  eps<-rnorm(n-1,0,sd)
  for(i in 2:n){
    X[i] <- rho*X[i-1]+eps[i-1]
  }
  return(X)
}



X <- simAR(10000, 0.7, sqrt(1-0.7^2)) #var = 0.09



(rho_1 <- est.rho.MLE(X))
(rho_2 <- est.rho.kim(X)$rho)
(sigma_1 <- 1-rho_1^2)
(sigma_2 <- est.rho.kim(X)$sigma)
