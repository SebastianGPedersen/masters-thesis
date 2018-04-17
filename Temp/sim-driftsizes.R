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
J <- sim.addjump(hest, alpha = 0.5, c_1 = 0.2, interval_length = 0.1)
J <- sim.path(path = 1, sim.data = J)

# bigfunction that returns both Raw/VB/VBDB
bursts <- function(settings, burstsetting, plot = F){
  alpha <- burstsetting$alpha
  beta <- burstsetting$beta
  burst_time <- burstsetting$burst_time
  interval_length <- burstsetting$interval_length
  c_1 <- burstsetting$c_1
  c_2 <- burstsetting$c_2
  
  sims<-sim.heston(setting)
  
  sims.vb<-sim.addvb(sims,    burst_time = burst_time, interval_length = interval_length,
                     c_2 = c_2, beta  = beta)
  
  sims.db<-sim.adddb(sims.vb, burst_time = burst_time, interval_length = interval_length,
                     c_1 = c_1,  alpha = alpha)
  
  
  #Get a single path
  path = 1
  
  Heston_path = sim.path(path,sims)$Y
  vb_path = sim.path(path, sims.vb)$Y
  vbdb_path = sim.path(path, sims.db)$Y
  
  if(!plot){
    return(list(raw = sim.path(path,sims), vb = sim.path(path, sims.vb), vbdb = sim.path(path, sims.db)))
  }
  
  #Set time index for plot
  time_index = sims$time/settings$mat
  
  
  #Only plot from 0.35 to 0.65
  index = ((time_index >=0.35) & (time_index <=0.65))
  Heston_path = Heston_path[index]
  vb_path = vb_path[index]
  vbdb_path = vbdb_path[index]
  time_index = time_index[index]
  
  df_y_values <- data.frame(cbind(Heston_path,vb_path, vbdb_path))
  names <- c("No burst","+Volatility burst", "+drift burst")
  
  tmp <- list()
  for (i in 1:3) {
    tmp[[i]] <- data.frame(
      x_akse = time_index,
      y_akse = df_y_values[,i],
      rx = names[i]
    )
  }
  tmp <- do.call(rbind, tmp)
  
  ggplot(tmp, aes(x_akse, y_akse,color = rx)) +
    geom_line() +
    scale_color_manual("", values = c("black", "red", "blue")) +
    xlab("time") + ylab("log-return")  + 
    ggtitle("Log-return of asset") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
}

burst<-sim.burstsetting(alpha = 0.6, beta = 0.2 ,c_1 = 0.2, c_2 = 0.03, interval_length = 0.1)

set.seed(1234)

sim.bursts <- bursts(setting, burst, F)

sims <- sim.bursts$raw # choose the one with vbdb

#sims <- J # USE JUMP INSTEAD


#plot(sims$Y, type = "l")

theta <- 1

k <- theta*152#1*sqrt(length(sims$time))
prev <- c(rep(0,k-2+2),est.NewPreAverage(diff(sims$Y), k)) # Kim does not divide by k - scaling does not change results

tind <- seq(from = 1000, to = 23000, by = 10)

# 1*52*7*24*60*60*1000 - YEAR TO MILLISECONDS
data <- list(time = sims$time*1*52*7*24*60*60*1000, Y = prev, raw = sims$Y)

dt <- diff(data$time)[1]
hd <- 300*dt
hv <- 1500*dt

conf = 0.95

# Estimation of mu/sig

mu<-est.mu.new(data = data, hd = hd, t.index = tind, kn = k)

sig <- est.sigma.new(data, hv=hv, t.index = tind, noisefun = est.noise.iid, theta = theta, kn = k)

# Calculate T
Tstat<-teststat(mu, sig, hd, hv)

# Calculate T*
Tstar<-tstar(Tstat)$tstar

# fit rho
rho <- est.rho(Tstat$test)

z<-est.z_quantile(rho$m, rho$rho, conf)$qZm
res<-Tstar>=z

a<-c(res, Tstar, z)
names(a) <- c("T/F", "Tstar", "z")

print(a)

plot(Tstat$test, type = "l")

# EXPORT TO MATLAB FOR PARALLEL VIEW
if(0 == 1){
  exp<-cbind(sims$time*1*52*7*24*60*60*1000,sims$Y)
  tim <- sims$time[tind]*1*52*7*24*60*60*1000
  
  write.csv(x = exp, file = "C:/Users/Frederik/Dropbox/Lspeciale/exp.csv")
  write.csv(x = tim, file = "C:/Users/Frederik/Dropbox/Lspeciale/tim.csv")
  write.csv(x = mu$mu, file = "C:/Users/Frederik/Dropbox/Lspeciale/mu.csv")
  write.csv(x = sig$sig, file = "C:/Users/Frederik/Dropbox/Lspeciale/sig.csv")
  write.csv(x = Tstat$test, file = "C:/Users/Frederik/Dropbox/Lspeciale/test.csv")
  #require(R.matlab)
  #writeMat(con="C:/Users/Frederik/Dropbox/Lspeciale/test.m", x=as.matrix(exp))
}
