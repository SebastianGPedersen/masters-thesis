
setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("kernels/kernels.R")


est.mu.legacy <- function(data, hd, kern, t.index){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # t.points might not work correctly....
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  # mode-handling
  t<-data$time[t.index]
  ind = t.index
  #t should now be data$time points
  
  #n = length(t)
  tt = length(t)
  n = length(data$time)
  mu = numeric(tt)
  
  dy = c(0, data$Y)
  for(j in 1:tt){
    mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[2:n])
  }
  # return list
  return(list(time = t, mu = mu))
}


hest <- sim.heston(sim.setup(Npath = 1))

data <- list(time = hest$time, Y = diff(hest$Y[1,]))
data2 <- list(time = hest$time, Y = diff(c(0,hest$Y[1,]))) # TO APPEND ZERO OR NOT TO
tind <- seq(1000, 23400, by = 5)
# TESTING GROUNDS #

#data <- list(time = 1:100*0.01, Y = 1:99)
#tind <- 99

mu1<-est.mu.legacy(data, 1, kern = kern.leftexp$kern, t.index = tind)$mu
mu2<-est.mu(data, hd = 1, t.index = tind)$mu
mu3<-est.mu.next(data, hd = 1, t.index = tind)$mu
mu4<-est.mu.next.cpp(data, hd = 1, t.index = tind)$mu

max(abs(mu1-mu2))
max(abs(mu2-mu3))
max(abs(mu3-mu4))

require(microbenchmark)

microbenchmark(est.mu.legacy(data, 1, kern = kern.leftexp$kern, t.index = tind),
               est.mu(data, hd = 1, t.index = tind),
               est.mu.next(data, hd = 1, t.index = tind),
               est.mu.next.cpp(data, hd = 1, t.index = tind),
               times = 1)

# SIGMA
sig1<-est.sigma(data, hv = 1, t.index = tind, lag = 10)$sig
sig2<-est.sigma.next(data, hv = 1, t.index = tind, lag = 10)$sig
sig3<-est.sigma.next.cpp(data, hv = 1, t.index = tind, lag = 10)$sig

max(abs(sig1-sig2))
max(abs(sig2-sig3))

microbenchmark(est.sigma(data, hv = 1, t.index = tind, lag = 10),
               est.sigma.next(data, hv = 1, t.index = tind, lag = 10),
               est.sigma.next.cpp(data, hv = 1, t.index = tind, lag = 10),
               times = 1)
