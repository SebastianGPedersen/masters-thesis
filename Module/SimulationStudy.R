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


setting<-sim.setup(Nsteps = 23400, Npath = 1) # use this such that we can add things to this "constructor" without changing code


# simulation
sim<-sim.heston(setting)

# pathwise
m<-dim(sim$Y)[1]

Tstar=numeric(m)

for(i in 1:m){
  # Extract
  simpath<-sim.path(i, sim)
  # Estimation of mu/sig
  mu<-est.mu(simpath, 120, kern.leftexp)
  sig <- est.sigma(simpath, 600, "auto", kern.leftexp, kern.parzen)
  
  # Calculate T
  Tstat<-teststat(mu, sig, 1200, kern.leftexp)
  
  # Calculate T*
  Tstar[i]<-tstar(Tstat,2)$tstar
}

plot(Tstar)

#require(profvis)
require(microbenchmark)

# 1 sti à 23 400 steps

sims<-sim.path(1,sim)
mu<-est.mu(sims, 120, kern.leftexp)
sig<-est.sigma(sims, 600, "auto", kern.leftexp, kern.parzen)

microbenchmark(sim.heston(setting), times = 10)                                               # 519 ms
microbenchmark(sim.path(1, sim), times = 10)                                                   
microbenchmark(est.mu(sims, 120, kern.leftexp), times = 1)                                    # 24 s
microbenchmark(est.sigma(sims, 600, "auto", kern.leftexp, kern.parzen), times = 1)            # 1900s (31min)
microbenchmark(teststat(mu, sig, 1200, kern.leftexp), times = 10)                             # 15 ms

