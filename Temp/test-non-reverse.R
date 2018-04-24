setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/pre-average.R")
source("kernels/kernels.R")
source("simulation/jumps.R")

setting <- sim.setup(Npath = 2, Nstep = 23400, omega = 0.0000225)
hest <- sim.heston(setting)

alpha = 0.7; c_1 = 0.1; burst_time = 0.5; interval_length = 0.1; beta = 0.4; c_2 = 0.1;

J <- sim.addjump(hest, alpha = alpha, c_1 = c_1, burst_time = burst_time, interval_length = interval_length)


sims.vb<-sim.addvb(hest,    burst_time = burst_time, interval_length = interval_length,
                   c_2 = c_2, beta  = beta, reverse = F)

sims.db<-sim.adddb(sims.vb, burst_time = burst_time, interval_length = interval_length,
                   c_1 = c_1,  alpha = alpha, reverse = F)

db <- sim.path(path = 1, sim.data = sims.db)

J <- sim.path(path = 1, sim.data = J)
plot(J$Y, type = "l", ylim = c(max(J$Y,db$Y),min(J$Y,db$Y)))
lines(db$Y, type = "l", col = "red")
