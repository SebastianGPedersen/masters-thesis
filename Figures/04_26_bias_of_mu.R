library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/BlackScholes.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 1.6*10^(-5) #What Mathias wrote
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
hd_list <- seq(1,15,by = 1)/ (60*24*7*52)
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
    mu_hat <- numeric(length = temp_paths)
    
    for (i in 1:temp_paths) {
      #i <- 1
      single_path <- list(Y = diff(BS$Y[i,]), time = BS$time)
      mu_hat[i] <- sqrt(hd_list[hd])/sqrt(K2*theta)*est.mu(data = single_path, hd = hd_list[hd], t.index = desired_index)$mu[1]
      
    }
    var_bias_temp[memory,hd] <- mean(mu_hat^2)
  }
}

#save mean(var)
var_bias <- numeric(length = length(hd_list))
noise <- var_bias

for (hd in 1:length(hd_list)) {
  var_bias[hd] <- mean(var_bias_temp[,hd])
  noise[hd] <- mean(1/hd_list[hd]*(1.6*10^(-5))^2) / sqrt(K2*theta^2)
}
var_bias
noise

data <- data.frame(hd = hd_list, var_bias = var_bias)

qplot(hd, var_bias, data = data, geom = "line")

##Calculate evt. mean(1/sqrt(h_n)*omega^2)
#################### PLOT ####################
