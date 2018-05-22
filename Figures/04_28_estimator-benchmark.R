# TESTING AND BENCHMARKING
setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/pre-average.R")
source("kernels/kernels.R")
source("simulation/jumps.R")

# SIMULATION 
setting <- sim.setup(Npath = 2, Nstep = 23400, omega = 0.0000225)
hest <- sim.heston(setting)
sim <- sim.path(path = 1, sim.data = hest)
sim <- list(time = sim$time, Y = diff(sim$Y), raw = diff(sim$Y))

# ESTIMATION SETTINGS 
n = 23400
mat <- 6.5/(24*7*52)#*52*7*24*60*60 #In years
dt <- mat/n #In years
hd <- 10^(-3)*dt^(1/3) #If years

# BENCHMARKERIA
require(ggplot2)
require(microbenchmark)

N<-23

time <- matrix(NA, nrow = N, ncol = 5)

index <- 1000*(1:N)
for(i in 1:N){
  tind <- 1000*i
  if(i == 1){
    prev.mu <- est.mu.next(sim, hd = hd, t.index = tind)
    prev.sig <- est.sigma.next(sim, hv = hd, t.index = tind, lag = 10)
    time[i,] <- c(tind, summary(microbenchmark(est.mu(sim, hd = hd, t.index = tind),
                                       est.mu.next(sim, hd = hd, t.index = tind),
                                       est.sigma(sim, hv = hd, t.index = tind, lag = 10),
                                       est.sigma.next(sim, hv = hd, t.index = tind, lag = 10),
                                       times = 10000, unit = "ms"))$mean)
  }
  else{
    time[i,] <- c(tind, summary(microbenchmark(est.mu(sim, hd = hd, t.index = tind),
                                       est.mu.next(sim, prevmu = prev.mu, hd = hd, t.index = tind),
                                       est.sigma(sim, hv = hd, t.index = tind, lag = 10),
                                       est.sigma.next(sim, prevsig = prev.sig, hv = hd, t.index = tind, lag = 10),
                                       times = 10000, unit = "ms"))$mean)
    prev.mu <- est.mu.next(sim, hd = hd, t.index = tind)
    prev.sig <- est.sigma.next(sim, hv = hd, t.index = tind, lag = 10)
  }
}
data<-list(index = time[, 1])

data$index <- time[,1]; data$mu <- time[,2]; data$mu.optim <- time[,3]; data$sig <- time[,4]; data$sig.optim <- time[,5]

data<-data.frame(data)

ggplot() +
  geom_line(data=data, aes(x=index, y=mu, col = "est.mu"), size = 1) +
  geom_line(data=data, aes(x=index, y=mu.optim, col = "est.mu.optim"), size = 1) +
  xlab("t.index") + ylab("Run time (ms)")  + 
  ggtitle("Computation time of mu-estimators")

ggplot() +
  geom_line(data=data, aes(x=index, y=sig, col = "est.sigma"), size = 1) +
  geom_line(data=data, aes(x=index, y=sig.optim, col = "est.sigma.optim"), size = 1) +
  xlab("t.index") + ylab("Run time (ms)")  + 
  ggtitle("Computation time of sigma-estimators")
