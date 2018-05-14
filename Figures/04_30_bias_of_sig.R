library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/BlackScholes.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")


#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 1.6*10^(-5)#*(100000/n) #What Mathias wrote
omega2 <- omega^2
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 100
sigma <- 0.0225

#Because of lack of memory, it is done in 10 loops
n_loops <- 1

#List to final values
hd_list <- seq(1,20,by = 1)/ (60*24*7*52)
sig_mean_temp <- matrix(nrow = n_loops,ncol = length(hd_list))
sig_var_temp <- sig_mean_temp


for (memory in 1:n_loops) {
  p0 <- Sys.time()
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = omega, Nsteps = n, Npath = temp_paths)
  
  #################### LOOP OVER H_D ####################
  
  #H_d from 1min to 15min
  
  for (hd in 1:length(hd_list)) {
    
    #hd <- 1
    desired_index <- n-1 #Takes last index, so K uses as many points as possible
    mu_hat <- numeric(length = temp_paths)
    sig_hat <- mu_hat
    T_hat <- mu_hat
    
    for (i in 1:temp_paths) {
      #i <- 1
      single_path <- list(Y = diff(BS$Y[i,]), time = BS$time)
      #sig_hat[i] <- 1/(K2*sigma^2) * est.sigma(single_path, hv = hd_list[hd], t.index = desired_index, lag = 20)$sig[1]
      sig_hat[i] <- 1/(K2*sigma^2) * est.sigma.next(data = single_path, hv = hd_list[hd], t.index = desired_index, lag = 20)$sig[1]
      
    }
    
    sig_mean_temp[memory,hd] <- mean(sig_hat)
    sig_var_temp[memory,hd] <- var(sig_hat)
  }
  print(Sys.time()-p0)
}

#save mean(var)
sig_mean <- numeric(length = length(hd_list))
sig_var <- sig_mean

for (hd in 1:length(hd_list)) {
  sig_mean[hd] <- mean(sig_mean_temp[,hd])
  sig_var[hd] <- mean(sig_var_temp[,hd])
}

data <- data.frame(hd = hd_list, 
                   target = (1:length(hd_list)*0)+1, 
                   mean = sig_mean, 
                   lower = sig_mean - sqrt(sig_var), 
                   upper = sig_mean + sqrt(sig_var))

ggplot() +
  geom_line(data=data, aes(x=hd, y=mean, col = "mean"), size = 1) +
  geom_line(data=data, aes(x=hd, y=lower, col = "lower"), size = 1) +
  geom_line(data=data, aes(x=hd, y=upper, col = "upper"), size = 1) +
  #geom_line(data=data, aes(x=hd, y=T_bias, col = "T_bias"), size = 1) +
  geom_line(data=data, aes(x=hd, y=target, col = "1"), size = 1) +
  xlab("Bandwidth, hd") + ylab("Normalized Variance")


##Calculate evt. mean(1/sqrt(h_n)*omega^2)
#################### PLOT ####################
