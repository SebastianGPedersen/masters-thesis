setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")

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

require(ggplot2)
setting <- sim.setup(Npath = 2, Nstep = 23400, omega = 0.0225/23400*10)

#burst<-sim.burstsetting(alpha = 0.5, beta = 0.2 ,c_1 = 0.2, c_2 = 0.03)

burst<-sim.burstsetting(alpha = 0.55, beta = 0.2 ,c_1 = 0.2, c_2 = 0.015)

#burst<-sim.burstsetting(alpha = 0.6, beta = 0.2 ,c_1 = 0.2, c_2 = 0.03)

set.seed(1234)

bursts(setting, burst, T)
sim.bursts <- bursts(setting, burst, F)

sims <- sim.bursts$vbdb

data <- est.EveryOtherDiffData(sims)

tind <- seq(from = 2000, to = 9000, by = 10)
#tind <- 23400/4

dt <- diff(sims$time)[1]
hd <- 100*dt
hv <- 1200*dt
conf = 0.95

# Estimation of mu/sig
# mu<-est.mu(data, hd, kern.leftexp, t.index = tind)
# sig <- est.sigma(data, hv, kern.leftexp, kern.parzen, t.index = tind, lag = "auto")

mu<-est.mu.next(data = data, hd = hd, t.index = tind)
sig <- est.sigma.next(data, hv = hv, t.index = tind, lag = "auto")

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
points(3500, Tstat$test[3500])
