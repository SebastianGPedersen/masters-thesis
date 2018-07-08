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
setting <- sim.setup(Npath = 249, Nsteps = 100000, omega = 1.6*10^-5)

seed<-2342
set.seed(seed)

# setup control and uneven
hest    <- sim.heston(setting)
hest.dt <- sim.heston.uneven(setting, seed = seed)

# LOG RETURNS
hest$Y <- t(diff(t(as.matrix(hest$Y))))
hest.dt$Y <- t(diff(t(as.matrix(hest.dt$Y))))

# ESTIMATION PARAMETERS
hd <- 300/(52*7*24*60*60) #(seconds)
hv <- 12*hd                                          
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)
t.index <- seq(offset, 100000, by = t.freq)

# ESTIMATION
# even
mu <- sqrt(hd)*est.mu.mat.2.0(data = hest, hd = hd, bandwidth_rescale = T)$mu[,t.index]#,t.index = t.index)$mu
sigma2 <- est.sigma.mat.2.0(data = hest, hv = hv, lag = lag, bandwidth_rescale = T)$sig[,t.index]#,t.index = t.index,lag = lag)$sig
Tstat <- mu/sqrt(sigma2)

# uneven.dt

mu.dt <- sqrt(hd)*est.mu.mat.2.0(data = hest.dt, hd = hd, bandwidth_rescale = T)$mu[,t.index]#,t.index = t.index)$mu
sigma2.dt <- est.sigma.mat.2.0(data = hest.dt, hv = hv, lag = lag, bandwidth_rescale = T)$sig[,t.index]#,t.index = t.index,lag = lag)$sig
Tstat.dt <- mu.dt/sqrt(sigma2.dt)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MEAN ACROSS DAYS
meanT<-colMeans(abs(Tstat))
meanT.dt<-colMeans(abs(Tstat.dt))

data<-data.frame(time = hest$time[t.index], Tstat = meanT, time.dt = hest.dt$time[t.index], Tstat.dt = meanT.dt)

require(ggplot2)
ggplot() +
  geom_point(data=data, aes(x=time*60*60, y=Tstat, col = "even dt"), size = 1) +
  geom_point(data=data, aes(x=time.dt*60*60, y=Tstat.dt, col = "uneven dt"), size = 1) +
  xlab("time") + ylab("Average absolute T value") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SIGMA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MEAN ACROSS DAYS
meanT<-colMeans(abs(sigma2))
meanT.dt<-colMeans(abs(sigma2.dt))

data<-data.frame(time = hest$time[t.index], Tstat = meanT, time.dt = hest.dt$time[t.index], Tstat.dt = meanT.dt)

require(ggplot2)
ggplot() +
  geom_point(data=data, aes(x=time*60*60, y=Tstat, col = "even dt"), size = 1) +
  geom_point(data=data, aes(x=time.dt*60*60, y=Tstat.dt, col = "uneven dt"), size = 1) +
  xlab("time") + ylab("Average absolute Sigma value") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MEAN ACROSS DAYS
meanT<-colMeans(abs(mu))
meanT.dt<-colMeans(abs(mu.dt))

data<-data.frame(time = hest$time[t.index], Tstat = meanT, time.dt = hest.dt$time[t.index], Tstat.dt = meanT.dt)

require(ggplot2)
ggplot() +
  geom_point(data=data, aes(x=time*60*60, y=Tstat, col = "even dt"), size = 1) +
  geom_point(data=data, aes(x=time.dt*60*60, y=Tstat.dt, col = "uneven dt"), size = 1) +
  xlab("time") + ylab("Average absolute Mu value") 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
