setwd(Sys.getenv("masters-thesis"))
source('Estimation/rho.R')
source("Estimation/gumbel.R")


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
T_star_func <- function(sim,m) {
  out <- max(abs(sim[(length(sim)-m):length(sim)]))
  return(out)
}



### Simulate T_stars
paths <- 500 #One second pr. 10 path
rho <- 0
m <- 2500
sd <- sqrt(1-rho^2)
a_m <- sqrt(2*log(m))
b_m <- a_m - 1/2*log(pi*log(m))/a_m

T_stars <- numeric(length = paths)
Gumbels <- T_stars
  
p0 <- Sys.time()
for (i in 1:paths) {
  sim <- simAR(10*m,rho,sd)
  T_stars[i] <- T_star_func(sim,m)
  Gumbels[i] <- (T_stars[i]-b_m)*a_m
}
print(Sys.time()-p0)

(my_T_quant <- quantile(T_stars,0.95)) #The simulated quantile of T -> 4.20
(my_gumbel_quant <- quantile(Gumbels,0.95)) #The simulated quantile of Gumpel -> 2.56
(theoretical_gumbel <- -log(-log(0.95))) #The theoretical level of Gumbel -> 2.97
(theoretical_T <- theoretical_gumbel/a_m+b_m) #The theoretical level of T^* -> 4.30


#Sandheden er at T aldrig kommer over 4.20


#### Prøv at lave reelle quantiles på max(T) ??
