#######################################################################
#                                                                     #
# FROM HERE WE CHAIN EVERYTHING TOGETHER TO FORM THE SIMULATION STUDY #
#                                                                     #
#######################################################################

# working directory should be /masters-thesis

setwd(Sys.getenv("masters-thesis"))

source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/teststat.R")
source("kernels/kernels.R")



study<-function(setting, t.index, h, q, beta, alpha){
  # simulation
  heston<-sim.heston(setting)
  
  # Add Burst
  Heston_vb <-   sim.addvb(Heston,    burst_time = 0.5, interval_length = 0.05, c_2 = 0.15, beta  = beta)
  Heston_vbdb <- sim.adddb(Heston_vb, burst_time = 0.5,   interval_length = 0.05, c_1 = 3,    alpha = alpha)
  
  Tstar = numeric(m)
  
  # pathwise
  m<-dim(sim$Y)[1]
  for(i in 1:m){
    # Extract
    simpath<-sim.path(i, sim)
    
    # Estimation of mu/sig
    mu<-est.mu(simpath, h, kern.leftexp, t.index = t.index)
    sig <- est.sigma(simpath, 5*h, kern.leftexp, kern.parzen, t.index = t.index, lag = "auto")
    
    # Calculate T
    Tstat<-teststat(mu, sig, h, kern.leftexp)
    
    # Calculate T*
    Tstar[i]<-tstar(Tstat)$tstar
    
    # plug in rest (q-quantile)
  }
  
  # return etwas
  return(list())
}

setting<-sim.setup(Nsteps = 23400, Npath = 1) # use this such that we can add things to this "constructor" without changing code

# find T points
tind<-seq(from = 60, to = 23400, by = 60)

# params

hset = c(120,300,600)/(3600*24*7*52)
alphaset = c(0.55, 0.65, 0.75)
betaset = c(0, 0.1, 0.2 ,0.3 ,0.4)

res = numeric(length(betaset)*length(alphaset)*length(hset)) # right?
for(beta in betaset){
  for(alpha in alphaset){
    for(h in hset)
    {
      study(setting, tind, h, alpha[1], beta[1])
    }
  }
}

sims<-sim.heston(sim.setup(Nsteps = 23400, Npath = 2))
sim<-sim.path(1, sims)

source("estimation/estimates.R")
mu<-est.mu(sim, 120/(3600*24*7*52), kern.leftexp, t.index = tind)
sigm<-est.sigma(sim, 5*120/(3600*24*7*52), kern.leftexp, kern.parzen, t.index = tind, lag = "auto")

plot(sigm$sig)
plot(test-sigm$sig)

plot(mu$mu, type = "l")
plot(sqrt(120/(3600*24*7*52))*mu$mu)
qqnorm(sqrt(120/(3600*24*7*52))*mu$mu)

