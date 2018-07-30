setwd(Sys.getenv("masters-thesis"))
source('Estimation/rho.R')

### Functions
simAR<-function(n, rho, sd){
  X<-rep(NA,n)
  X[1]<-rnorm(1, 0, sd/(1-rho))
  eps<-rnorm(n-1,0,sd)
  for(i in 2:n){
    X[i] <- rho*X[i-1]+eps[i-1]
  }
  return(X)
}


### Simulate T'es
paths <- 500 #One second pr. 10 path
rho <- 0.7
m <- 2500
sd <- sqrt(1-rho^2)
a_m <- sqrt(2*log(m))
b_m <- a_m - 1/2*log(pi*log(m))/a_m


rhos_kim <- numeric(length = paths)
rhos_mle <- rhos_kim
rhos_third_degree <- rhos_kim

p0 <- Sys.time()
for (i in 1:paths) {
  sim <- simAR(m,rho,sd)
  rhos_kim[i] <- est.rho.kim(sim)
  rhos_mle[i] <- est.rho.MLE(sim)
  rhos_third_degree[i] <- est.rho.third_degree(sim)
}

### mean
print(mean(rhos_kim[i]))
print(mean(rhos_mle[i]))
print(mean(rhos_third_degree[i])) #helt helt hen i vejret lige nu

### max abs distance
print(max(abs(rhos_kim-rhos_mle)))
print(max(abs(rhos_kim-rhos_third_degree)))
print(max(abs(rhos_mle-rhos_third_degree)))

