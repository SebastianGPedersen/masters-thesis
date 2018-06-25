setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/pre-average.R")
source("kernels/kernels.R")

# SIMULATE
setting <- sim.setup(Npath = 2, Nsteps = 50000, omega = 1.6*10^-5)
seed<-2342
set.seed(seed)
hest    <- sim.heston(setting)
hest.dt <- sim.heston.uneven(setting, seed = seed)

# ADD burst
Heston_vb   <- sim.addvb(hest, burst_time = 0.5, interval_length = 0.05, c_2 =) #WHAT IS C2?, beta = 0.1, reverse = F)
Heston_vbdb <- sim.adddb(Heston_vb, burst_time=0.5, interval_length=0.05, c_1 = , alpha=0.65, reverse = F)

Heston_vb.dt   <- sim.addvb(hest.dt, burst_time = 0.5, interval_length = 0.05, c_2 = , beta = 0.1, reverse = F)
Heston_vbdb.dt <- sim.adddb(Heston_vb.dt, burst_time=0.5, interval_length=0.05, c_1 = , alpha=0.8, reverse = F)

hest <- Heston_vbdb
hest.dt <- Heston_vbdb.dt


data<-list(time = hest$time, Y = hest$Y[1,], time2 = hest.dt$time, Y2 = hest.dt$Y[1,])

data<-data.frame(data)

require(ggplot2)
ggplot() +
  geom_line(data=data, aes(x=time, y=Y, col = "even dt"), size = 1) +
  geom_line(data=data, aes(x=time2, y=Y2, col = "uneven dt"), size = 1) +
  xlab("time") + ylab("log-return") 

# ESTIMATION
n = 50000
mat <- 6.5/(24*7*52)#*52*7*24*60*60 #In years
dt <- mat/n #In years
hd <- 10^(-3)*dt^(1/3) #If years
tind <- seq(10000, 40000, by = 10)

hest<-list(time = hest$time, Y = diff(hest$Y[1,]))
hest.dt<-list(time = hest.dt$time, Y = diff(hest.dt$Y[1,]))

# Plot T
mu <- est.mu.next(data = hest, hd = hd, t.index = tind)
mu2 <- est.mu.next(data = hest.dt, hd = hd, t.index = tind)
K2 <- 0.5
test <- list(time = mu$time, test = sqrt(hd/K2)*mu$mu/sqrt(est.sigma.raw.next(data = hest, hv = hd, t.index = tind)$sig))
test.dt <- list(time = mu2$time, test = sqrt(hd/K2)*mu2$mu/sqrt(est.sigma.raw.next(data = hest.dt, hv = hd, t.index = tind)$sig))

plot.data<-data.frame(test)
plot.data.dt<-data.frame(test.dt)

require(ggplot2)
if(F){
  ggplot(plot_data_frame) + 
    stat_qq(aes(sample = T_estimator, colour = Bandwidth)) +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Quantile of standard normal") + ylab("Quantile of T")
}
