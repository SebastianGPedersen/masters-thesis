library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/BlackScholes.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 1.6*10^(-5)*(100000/n) #What Mathias wrote
omega2 <- omega^2
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 1000
theta <- 0.0225

#Because of lack of memory, it is done in 10 loops
n_loops <- 10

#List to final values
hd_list <- seq(1,20,by = 1)/ (60*24*7*52)
var_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))

for (memory in 1:n_loops) {
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  BS <- sim.BlackScholes(mean = 0, sd = theta, omega = omega, Nsteps = n, Npath = temp_paths)
  
  #################### LOOP OVER H_D ####################
  
  #H_d from 1min to 15min
  
  for (hd in 1:length(hd_list)) {
    
    #hd <- 1
    desired_index <- n-1 #Takes last index, so K uses as many points as possible
    sig_hat <- numeric(length = temp_paths)
    
    for (i in 1:temp_paths) {
      #i <- 1
      single_path <- list(Y = diff(BS$Y[i,]), time = BS$time)
      sig_hat <- 1/theta^2 * est.sigma(single_path, hv = hd_list[hd], t.index = desired_index, lag = 10)$sig
    }
    var_bias_temp[memory,hd] <- mean(sig_hat)
  }
}

#save mean(var)
var_bias <- numeric(length = length(hd_list))
noise <- var_bias

for (hd in 1:length(hd_list)) {
  var_bias[hd] <- mean(var_bias_temp[,hd])
  noise[hd] <- mean(1/hd_list[hd]*omega^2) / sqrt(K2*theta^2)
}

data <- data.frame(hd = hd_list, target = (1:length(hd_list)*0)+1, var_bias = var_bias)

ggplot() +
  geom_line(data=data, aes(x=hd, y=var_bias, col = "est.sig/sig"), size = 1) +
  geom_line(data=data, aes(x=hd, y=target, col = "target"), size = 1) +
  xlab("Bandwidth, hd") + ylab("Normalized Variance")


##Calculate evt. mean(1/sqrt(h_n)*omega^2)
#################### PLOT ####################
