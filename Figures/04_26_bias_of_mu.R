library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 0 #What Mathias wrote
omega2 <- omega^2
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 1000

#List to final values
hd_list <- seq(1,15,by = 1)/ (60*24*7*52)
var_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))


#Because of lack of memory, it is done in 10 loops
n_loops <- 10

for (memory in 1:n_loops) {
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega, xi = 0) #6.5 hours
  Heston <- sim.heston(settings)
  
  #################### LOOP OVER H_D ####################
  
  #H_d from 1min to 15min
    
  for (hd in 1:length(hd_list)) {
    
    #hd <- 1
    desired_index <- n #Takes last index, so K uses as many points as possible
    mu_hat <- numeric(length = temp_paths)
    
    for (i in 1:temp_paths) {
      #i <- 1
      single_path <- list(Y = Heston$Y[i,], time = Heston$time, raw = Heston$raw[i,])
      mu_hat[i] <- sqrt(hd_list[hd])*est.mu(data = single_path, hd,t.index = desired_index)$mu[1]
      mu_hat[i] <- mu_hat[i] / sqrt(K2*settings$theta) #Scale to supposably N(0,1). Theta = sigma^2
    }
    
    var_bias_temp[memory,hd] <- mean(mu_hat^2)
  }
}


#save mean(var)
var_bias <- numeric(length = length(hd_list))
noise <- var_bias

for (hd in 1:length(hd_list)) {
  var_bias[hd] <- mean(var_bias_temp[,hd])
  noise[hd] <- mean(1/hd_list[hd]*(1.6*10^(-5))^2) / sqrt(K2*settings$theta^2)
}
var_bias
noise

##Calculate evt. mean(1/sqrt(h_n)*omega^2)
#################### PLOT ####################
