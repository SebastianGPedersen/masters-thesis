library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/BlackScholes.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

p0 <- Sys.time()

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 1.6*10^(-5)#*(100000/n) #What Mathias wrote
omega2 <- omega^2
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 500
sigma <- 0.0225
(noise_ratio <- omega/sigma*sqrt(n))

#Alternative
#noise_ratio <- 1
#omega <- noise_ratio*sigma/sqrt(n)

#Lag length choice with (3) in Barndorff-Nielsen (2009)
c_star <- 3.5134
xi <- sqrt(omega^2/sqrt(mat^2*sigma^4))
lag <- c_star*xi^(4/5)*n^(3/5)

#Because of lack of memory, it is done in 10 loops
n_loops <- 5

#List to final values
h_mu <- 5 / (60*24*7*52)
ratio_list <- seq(2,20,by = 1) / 2
T_above_heston <- matrix(nrow = n_loops,ncol = length(ratio_list))
T_above_small_jump <- T_above_heston
T_above_large_jump <- T_above_heston
var_T_bias_temp <- T_above_heston

for (memory in 1:n_loops) {
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  
  #Heston simulations
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  Heston_jump_small <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = 0.3, alpha = 0.55)
  Heston_jump_large <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = 0.016, alpha = 0.8)
  
  #BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = omega, Nsteps = n, Npath = temp_paths)
  
  #################### LOOP OVER H_D ####################
  
  #H_d from 1min to 15min
  
  for (hd in 1:length(ratio_list)) {
    
    #hd <- 1
    desired_index <- n/2 #
    mu_hat <- numeric(length = temp_paths)
    sig_hat <- mu_hat
    T_hat_heston <- mu_hat
    T_hat_small_jump <- mu_hat
    T_hat_large_jump <- mu_hat
    
    for (i in 1:temp_paths) {
      #i <- 1
      
      #Heston
      single_path <- list(Y = diff(Heston$Y[i,]), time = Heston$time)
      sig_hat[i] <- est.sigma(single_path, hv = h_mu*ratio_list[hd], t.index = desired_index, lag = lag)$sig[1]
      mu_hat[i] <- est.mu(data = single_path, hd = h_mu, t.index = desired_index)$mu[1]
      T_hat_heston[i] <- sqrt(h_mu)*mu_hat[i]/sqrt(sig_hat[i]) #K2*sigma cancels out and we end with sqrt(h_n)*\hat{\mu} / sqrt(\hat{\Sigma}^2)
      
      #Small jump
      single_path <- list(Y = diff(Heston_jump_small$Y[i,]), time = Heston$time)
      sig_hat[i] <- est.sigma(single_path, hv = h_mu*ratio_list[hd], t.index = desired_index, lag = lag)$sig[1]
      mu_hat[i] <- est.mu(data = single_path, hd = h_mu, t.index = desired_index)$mu[1]
      T_hat_small_jump[i] <- sqrt(h_mu)*mu_hat[i]/sqrt(sig_hat[i]) #K2*sigma cancels out and we end with sqrt(h_n)*\hat{\mu} / sqrt(\hat{\Sigma}^2)
      
      #Large_jump
      single_path <- list(Y = diff(Heston_jump_large$Y[i,]), time = Heston$time)
      sig_hat[i] <- est.sigma(single_path, hv = h_mu*ratio_list[hd], t.index = desired_index, lag = lag)$sig[1]
      mu_hat[i] <- est.mu(data = single_path, hd = h_mu, t.index = desired_index)$mu[1]
      T_hat_large_jump[i] <- sqrt(h_mu)*mu_hat[i]/sqrt(sig_hat[i]) #K2*sigma cancels out and we end with sqrt(h_n)*\hat{\mu} / sqrt(\hat{\Sigma}^2)
    }
    
    #Safe perc above 0.95
    interval_196 <- qnorm(0.975)
    T_above_heston[memory,hd] <- sum((abs(T_hat_heston) > interval_196))/length(T_hat_heston)
    T_above_small_jump[memory,hd] <- sum((abs(T_hat_small_jump) > interval_196))/length(T_hat_small_jump)
    T_above_large_jump[memory,hd] <- sum((abs(T_hat_large_jump) > interval_196))/length(T_hat_large_jump)
    
    var_T_bias_temp[memory,hd] <- mean(T_hat^2)
  }
}

#Save mean(T_above_percentage)
Heston_above <- numeric(length = length(ratio_list))
small_above <- Heston_above
large_above <- Heston_above
var_T_bias <- Heston_above

for (hd in 1:length(ratio_list)) {
  Heston_above[hd] <- mean(T_above_heston[,hd])
  small_above[hd] <- mean(T_above_small_jump[,hd])
  large_above[hd] <- mean(T_above_large_jump[,hd])
  var_T_bias[hd] <- mean(var_T_bias_temp[,hd])
}

print(var_T_bias) #To see if unit variance 

#################### PLOT THE PERCENTAGE ABOVE ####################
data <- data.frame(hd = ratio_list, 
                   target = (1:length(ratio_list)*0)+0.05, 
                   Heston_above = Heston_above,
                   small_above = small_above,
                   large_above = large_above)

ggplot() +
  geom_line(data=data, aes(x=hd, y=target, col = "5 percent"), size = 1) +
  geom_line(data=data, aes(x=hd, y=Heston_above, col = "Heston"), size = 1) +
  geom_line(data=data, aes(x=hd, y=small_above, col = "Small jump"), size = 1) +
  geom_line(data=data, aes(x=hd, y=large_above, col = "Large jump"), size = 1) +
  xlab("Ratio of h_sigma / h_mu") + ylab("Percentage above T_95")

print(Sys.time()-p0)

