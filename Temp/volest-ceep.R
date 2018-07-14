setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/BlackScholes.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("kernels/kernels.R")
source("vol/vol_estimators.R")


# RCPP TEST ON ESTIMATES #
require(microbenchmark)
require(Rcpp)

sourceCpp("temp/volest.cpp")

# settings
K <- function(theta, N){
  return(floor(theta*sqrt(N))+floor(theta*sqrt(N))%%2)
}

Nsteps <- 10000

# simulation
hest <- sim.heston(sim.setup(Nsteps = Nsteps, Npath = 2))

Y <- diff(hest$Y[1,])

pa <- vol.est.preA(Y, K = K(1, Nsteps+1))

omega2 <- vol.est.omega2(Y)

RV <- vol.est.RVstar(pa, K(1, Nsteps+1), omega2Est = omega2, theta = 1)

BV <- vol.est.BVstar(pa, K(1, Nsteps+1), omega2Est = omega2, theta = 1)

#pa-vol_est_preA_cpp(Y, K(1, Nsteps+1)) #checkeria

RV-vol_est_RVstar(pa, K(1, Nsteps+1), omega2, theta = 1)

BV-vol_est_BVstar(pa, K(1, Nsteps+1), omega2, theta = 1)

vol.est.BV(pa)-vol_est_BV(pa)

microbenchmark(vol.est.preA(Y, K = K(1, Nsteps+1)),
               vol_est_preA_cpp(Y, K(1, Nsteps+1)),
               vol.est.RVstar(pa, K(1, Nsteps+1), omega2Est = omega2, theta = 1),
               vol_est_RVstar(pa, K(1, Nsteps+1), omega2, theta = 1),
               vol.est.BVstar(pa, K(1, Nsteps+1), omega2Est = omega2, theta = 1),
               vol_est_BVstar(pa, K(1, Nsteps+1), omega2, theta = 1))



