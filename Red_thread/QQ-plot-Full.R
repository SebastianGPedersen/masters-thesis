setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/rescaling.R")
source("kernels/kernels.R")
source("spy/dataFunctions.R")
source("spy/datahandling.R")

library(ggplot2)
source("Vol/BSS_Sim.R")
source("Simulation/Noise.R")
source("Estimation/estimates_reloaded.R")

# IMPORT BSS as "HESTON"
setwd(Sys.getenv("masters-thesis-data"))
load("BSSsim.Rdata") # called BSSsim
setwd(Sys.getenv("masters-thesis"))

setting <- sim.setup(Nstep = 20000, Npath = 1000)

BSSsim$Y <- BSSsim$X + setting$omega*rnorm(   length(BSSsim$time), 0 , 1   ) # add microstruct

# LOG RETURNS
BSSsim$Y <- t(diff(t(as.matrix(BSSsim$Y))))


# ESTIMATION PARAMETERS
hd <- 300/(52*7*24*60*60) #(seconds)
hv <- 12*hd                                          
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)
t.index <- seq(offset, setting$Nsteps, by = t.freq)

# ESTIMATION
# even
mu <- sqrt(hd)*est.mu.mat.2.0(data = BSSsim, hd = hd, bandwidth_rescale = T)$mu[,t.index]#,t.index = t.index)$mu
sigma2 <- est.sigma.mat.2.0(data = BSSsim, hv = hv, lag = lag, bandwidth_rescale = T)$sig[,t.index]#,t.index = t.index,lag = lag)$sig
Tstat <- mu/sqrt(sigma2)

plot_data_frame <- data.frame(list(T_estimator = as.numeric(Tstat)))

ggplot(plot_data_frame) + 
  stat_qq(aes(sample = T_estimator)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Quantile of standard normal") + ylab("Quantile of T") +
  theme(legend.position=c(0.87,0.16))


