set.seed(100)
library(ggplot2)
library(grid)
library(gridExtra)

setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/BlackScholes.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")
source("Estimation/estimates_reloaded.R")
source("Estimation/estimates_revolution.R")

p0 <- Sys.time()

#################### PARAMETERS THAT DON'T CHANGE ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400 /7 #Cheating again
mat <- 6.5/(24*7*52) /7
dt <- mat/n
Npaths <- 5000
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66
lag <- 10
ratio <- 15

#Because of lack of memory, it is done in loops
n_loops <- Npaths/100
desired_index <- n-1

#List to final values
hd_list <- c(2,5,10) / (60*24*7*52) #2min, 5min and 10min
T_estimator <- matrix(nrow = Npaths,ncol = length(hd_list))

for (memory in 1:n_loops) {
  #memory <- 1
  print(paste("memory",memory,"out of",n_loops))
  temp_paths <- Npaths / n_loops
  
  #Heston simulations
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  
  #### Calculate dY ####
  dY <- t(diff(t(Heston$Y))) #Has to be transposed for 'diff' to work as desired
  data <- list(Y = dY, time = Heston$time)    
  
  for (hd in 1:length(hd_list)) {
    sig_hat <- est.sigma.mat.2.0(data, hv = ratio*hd_list[hd], lag = lag)$sig[,desired_index]
    mu_hat <- est.mu.mat.2.0(data, hd = hd_list[hd])$mu[,desired_index]
    T_estimator[((memory-1)*temp_paths+1):(memory*temp_paths),hd] <- sqrt(hd_list[hd]) * mu_hat / sqrt(sig_hat)    
  }
}


#################### PLOT ####################
#Re-shape to data-frame
plot_data_frame <- data.frame(do.call("rbind",
                           list(cbind(T_estimator[,1],rep(" 2min",dim(T_estimator)[1])), 
                                cbind(T_estimator[,2],rep(" 5min",dim(T_estimator)[1])),
                                cbind(T_estimator[,3],rep("10min",dim(T_estimator)[1])))))

colnames(plot_data_frame) <- c("T_estimator", "Bandwidth")

#Numeric and factor
plot_data_frame$T_estimator <- as.numeric(as.character(plot_data_frame$T_estimator))
plot_data_frame$Bandwidth <- as.factor(plot_data_frame$Bandwidth)


g2 <- ggplot(plot_data_frame) + 
  stat_qq(aes(sample = T_estimator, colour = Bandwidth)) +
  geom_abline(intercept = 0, slope = 1) +
  #ylim(-5,5) +
  xlim(-4,4) +
  xlab("Quantile of standard normal") + ylab("Quantile of T") +
  theme(legend.position=c(0.87,0.16))
g2


