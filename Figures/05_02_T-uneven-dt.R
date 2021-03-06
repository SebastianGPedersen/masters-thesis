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
setting <- sim.setup(Npath = 2, Nsteps = 23399, omega = 0)

hest    <- sim.heston(setting)
hest.dt <- sim.heston.uneven(setting)

# --- LOOK AT TRADES DIST - REAL/HERE--- #
# REAL DATA
load("Simulation/trades.RDa")
data <- data.frame(Time = trades$time, Trades = trades$trades)

# OURS
s_dt <-hest.dt$time*52*7*24*60*60
s <- hest$time*52*7*24*60*60


time<-as.POSIXct.numeric(s, tz = "UTC", origin = "2014-01-02 09:30:00 UTC")
buckets<-seq(time[1], time[length(time)], by = 300  )
breaks <- .bincode(time, breaks = c(buckets, time[length(time)]) )
sim<-data.table(Time = time, group = c(1,breaks[2:length(breaks)])) # prevent boundary fuckery
N <- sim[, .N, by = group]$N


time_dt<-as.POSIXct.numeric(s_dt, tz = "UTC", origin = "2014-01-02 09:30:00 UTC")
#buckets_dt<-seq(time_dt[1], time_dt[length(time_dt)], by = 300  )
breaks_dt <- .bincode(time_dt, breaks = buckets)
dt<-data.table(Time = time_dt, group = c(1,breaks_dt[2:length(breaks_dt)])) # prevent boundary fuckery
N_dt <- dt[, .N, by = group]$N

plot.data <- data.frame(Time = buckets[-1], N_dt = N_dt, N = N)

# ----------------------- QQ PLOT SECTION -------------------------

require(ggplot2)
require(grid)
require(gridExtra)

p0 <- Sys.time()

omega2 <- 2.64*10^(-10)*25
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400 /7
mat <- 6.5/(24*7*52) /7
dt <- mat/n
Npaths <- 1000
sigma2 <- 0.0457/25
sigma <- sqrt(sigma2)
bandwidth_ratio <- 1
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66
lag <- 100

#Because of lack of memory, it is done in loops
n_loops <- 10

#List to final values
h_list <- c(2,5,10) / (60*24*7*52) #2min, 5min and 10min
T_estimator <- matrix(nrow = Npaths,ncol = length(h_list))


for (memory in 1:n_loops) {
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  
  #Heston simulations
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston.uneven(settings)
  #Heston <- sim.heston(settings)
  
  for (h_mu in 1:length(h_list)) {
    #h_mu <- 1
    desired_index <- n-1 #Just takes last index so K has as many obs as possible
    
    for (i in 1:temp_paths) {
      #i <- 1
      #Heston
      single_path <- list(Y = diff(Heston$Y[i,]), time = Heston$time)
      sig_hat <- est.sigma(single_path, hv = h_list[h_mu]*bandwidth_ratio, t.index = desired_index, lag = lag)$sig[1]
      mu_hat <- est.mu(data = single_path, hd = h_list[h_mu], t.index = desired_index)$mu[1]
      T_estimator[(memory-1)*temp_paths+i,h_mu] <- sqrt(h_list[h_mu])*mu_hat/sqrt(sig_hat) #K2*sigma cancels out and we end with sqrt(h_n)*\hat{\mu} / sqrt(\hat{\Sigma}^2)
    }
  }
}


#Re-shape to data-frame
plot_data_frame <- data.frame(do.call("rbind",
                                      list(cbind(T_estimator[,1],rep(" 2min",dim(T_estimator)[1])), 
                                           cbind(T_estimator[,2],rep(" 5min",dim(T_estimator)[1])),
                                           cbind(T_estimator[,3],rep("10min",dim(T_estimator)[1])))))

colnames(plot_data_frame) <- c("T_estimator", "Bandwidth")

#Numeric and factor
plot_data_frame$T_estimator <- as.numeric(as.character(plot_data_frame$T_estimator))
plot_data_frame$Bandwidth <- as.factor(plot_data_frame$Bandwidth)

# ------

# TRIM TIME
timebreaks <- c(buckets[1], buckets[floor(length(buckets)/2)], buckets[length(buckets)])

# PLOT
require(ggplot2)
g1 <- ggplot(data = plot.data) +
  geom_point(aes(x = Time, y = N, group = 1, col = "even")) +
  geom_point(aes(x = Time, y = N_dt, group = 1, col = "uneven")) +
  ylab("Trades per 5min bucket") +
  scale_x_datetime("Time bucket", breaks = timebreaks)
g1

g2 <- ggplot(plot_data_frame) + 
  stat_qq(aes(sample = T_estimator, colour = Bandwidth)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Quantile of standard normal") + ylab("Quantile of T")
g2
grid.arrange(g1,g2,nrow = 1,
             top = textGrob("",gp=gpar(fontsize=20)))
