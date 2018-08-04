setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
source("estimation/estimates_reloaded.R")
source("estimation/estimates_revolution.R")

# CHECK FOR SEASONALITY OF T IN UN-EVEN DT STUDY TO DETERMINE ORIGIN OF SEASONALITY

# SIMULATE
setting <- sim.setup(Npath = 1000, Nsteps = 20000, omega = 1.6*10^-5)

seed<-2342
set.seed(seed)

# setup control and uneven
hest.dt <- sim.heston.uneven(setting, seed = seed)
# neven

# IMPORT BSS as "HESTON"
#setwd(Sys.getenv("masters-thesis-data"))
load("BSSsim_uneven.Rdata") # called BSSsim
setwd(Sys.getenv("masters-thesis"))

BSSsim_uneven$Y <- BSSsim_uneven$X + setting$omega*rnorm(   length(BSSsim_uneven$time), 0 , 1   ) # add microstruct

# LOG RETURNS
BSSsim_uneven$Y <- t(diff(t(as.matrix(BSSsim_uneven$Y))))
hest.dt$Y <- t(diff(t(as.matrix(hest.dt$Y))))

# ESTIMATION PARAMETERS
hd <- 300/(52*7*24*60*60) #(seconds)
hv <- 12*hd                                          
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)
t.index <- seq(offset, setting$Nsteps, by = t.freq)

# ESTIMATION
# even
mu <- sqrt(hd)*est.mu.mat.2.0(data = BSSsim_uneven, hd = hd, bandwidth_rescale = T)$mu[,t.index]#,t.index = t.index)$mu
sigma2 <- est.sigma.mat.2.0(data = BSSsim_uneven, hv = hv, lag = lag, bandwidth_rescale = T)$sig[,t.index]#,t.index = t.index,lag = lag)$sig
Tstat <- mu/sqrt(sigma2)

# uneven.dt

mu.dt <- sqrt(hd)*est.mu.mat.2.0(data = hest.dt, hd = hd, bandwidth_rescale = T)$mu[,t.index]#,t.index = t.index)$mu
sigma2.dt <- est.sigma.mat.2.0(data = hest.dt, hv = hv, lag = lag, bandwidth_rescale = T)$sig[,t.index]#,t.index = t.index,lag = lag)$sig
Tstat.dt <- mu.dt/sqrt(sigma2.dt)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MEAN ACROSS DAYS
meanT<-apply(Tstat, MARGIN = 2, FUN = var)
meanT.dt<-apply(Tstat.dt, MARGIN = 2, FUN = var)

data<-data.frame(time = BSSsim_uneven$time[t.index], Tstat = meanT, time.dt = hest.dt$time[t.index], Tstat.dt = meanT.dt)

require(ggplot2)
ggplot() +
  geom_point(data=data, aes(x=time*60*60, y=Tstat, col = "Non-equidistant BSS"), size = 1) +
  geom_point(data=data, aes(x=time.dt*60*60, y=Tstat.dt, col = "Non-equidistant Heston"), size = 1) +
  geom_hline(yintercept = 1) +
  xlab("Time") + ylab("Average variance of the  T-estimator") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SIGMA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MEAN ACROSS DAYS
meanT<-colMeans(sigma2)
meanT.dt<-colMeans(sigma2)

data<-data.frame(time = BSSsim_uneven$time[t.index], Tstat = meanT, time.dt = hest.dt$time[t.index], Tstat.dt = meanT.dt)

require(ggplot2)
ggplot() +
  geom_point(data=data, aes(x=time*60*60, y=Tstat, col = "Non-equidistant BSS"), size = 1) +
  geom_point(data=data, aes(x=time.dt*60*60, y=Tstat.dt, col = "Non-equidistant Heston"), size = 1) +
  xlab("Time") + ylab("Average Sigma estimator") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MUSIG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MEAN ACROSS DAYS
meanmu<-apply(mu, MARGIN = 2, FUN = var)
meanmu.dt<-apply(mu.dt, MARGIN = 2, FUN = var)

data<-data.frame(time = BSSsim_uneven$time[t.index], Mu = meanmu,
                 time.dt = hest.dt$time[t.index], Mu.dt = meanmu.dt,
                 Sig = meanT, Sig.dt = meanT.dt)

require(ggplot2)
ggplot() +
  geom_point(data=data, aes(x=time*60*60, y=Mu, col = "Variance of drift estimator"), size = 1) +
  geom_point(data=data, aes(x=time*60*60, y=Sig, col = "Average volatility estimate"), size = 1) +
  xlab("time") + ylab("Value of estimator for BSS") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
