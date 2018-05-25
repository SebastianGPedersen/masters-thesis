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
setting <- sim.setup(Npath = 2, Nstep = 100, omega = 0.0000225)
hest <- sim.heston(setting)

data <- list(time = hest$time, Y = hest$Y[1,])
if(F){
  #mu
  est.mu.mat(hest, hd = 0.001, t.index = 50:55)
  est.mu(data, hd = 0.001, t.index = 50:55)
  est.mu.next(data, hd = 0.001, t.index = 50:55)
  
  # sigma
  est.sigma.mat(hest, hv = 0.001, t.index = 50:55, lag = 10)$sig
  est.sigma(data, hv = 0.001, t.index = 50:55, lag = 10)$sig
  est.sigma.next(data, hv = 0.001, t.index = 50:55, lag = 10)$sig
}

#mu
est.mu.mat(hest, hd = 0.001, t.index = 50:55)
est.mu.mat.next(hest, hd = 0.001, t.index = 50:55)

# sigma
est.sigma.mat(hest, hv = 0.001, t.index = 50:55, lag = 10)
est.sigma.mat.next(hest, hv = 0.001, t.index = 50:55, lag = 10)
