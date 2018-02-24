#######################################################################
#                                                                     #
# FROM HERE WE CHAIN EVERYTHING TOGETHER TO FORM THE SIMULATION STUDY #
#                                                                     #
#######################################################################

# working directory should be /masters-thesis

setwd("C:/Users/Frederik/Dropbox/Lspeciale/masters-thesis") # generalize

source("simulation/heston.R")
source("estimation/estimates.R")
source("estimation/teststat.R")
source("kernels/kernels.R")

#require(profvis)
require(microbenchmark)

setting<-sim.setup() # use this such that we can add things to this "constructor" without changing code

# simulation
sim<-sim.heston(setting)

# pathwise
m<-1#dim(sim$Y)[1]

for(i in 1:m){
  # Extract
  simpath<-sim.path(i, sim)
  # Estimation of mu/sig
  mu<-est.mu(simpath, 120, kern.leftexp)
  sig <- est.sigma(simpath, 600, "auto", kern.leftexp, kern.parzen)
  
  # Calculate T
  Tstat<-teststat(mu, sig, 1200, kern.leftexp)
  
  # Calculate T*
  Tstar<-tstar(Tstat,2)
}

microbenchmark(sim.path(1, sim))
microbenchmark(est.mu(simpath, 120, kern.leftexp))  # 50 ms
microbenchmark(est.sigma(simpath, 600, "auto", kern.leftexp, kern.parzen)) #4s with 10 lag and 1001 observations fuck
microbenchmark(teststat(mu, sig, 1200, kern.leftexp)) #185 µs

