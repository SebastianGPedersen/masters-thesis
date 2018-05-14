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
omega <- 1.6*10^(-5)*(100000/n) #What Mathias wrote
omega2 <- omega^2
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 200
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
n_loops <- 4

#List to final values
hd_list <- seq(1,20,by = 1)/ (60*24*7*52)
var_mu_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))
sig_temp <- matrix(nrow = n_loops,ncol = length(hd_list))
var_T_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))


for (memory in 1:n_loops) {
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
      sig_hat[i] <- 1/(K2*sigma^2) * est.sigma(single_path, hv = hd_list[hd], t.index = desired_index, lag = lag)$sig[1]
      mu_hat[i] <- sqrt(hd_list[hd])/sqrt(K2*sigma^2)*est.mu(data = single_path, hd = hd_list[hd], t.index = desired_index)$mu[1]
      T_hat[i] <- mu_hat[i]/sqrt(sig_hat[i]) #K2*sigma cancels out and we end with sqrt(h_n)*\hat{\mu} / sqrt(\hat{\Sigma}^2)
    }
    
    var_mu_bias_temp[memory,hd] <- mean(mu_hat^2) #Vi ved mean er 0, så korrekt Var er mean^2
    sig_temp[memory,hd] <- mean(sig_hat)
    print(length(T_hat)-length(na.omit(T_hat))) #Check how many is missing
    var_T_bias_temp[memory,hd] <- var(na.omit(T_hat)) #Mean er ikke 0, så korrekt var er 
  }
}

#################### TAKE MEAN OVER LOOPS  ####################
var_mu_bias <- numeric(length = length(hd_list))
sig_bias <- var_mu_bias
var_T_bias <- var_mu_bias
noise <- var_mu_bias

for (hd in 1:length(hd_list)) {
  var_mu_bias[hd] <- mean(var_mu_bias_temp[,hd])
  sig_bias[hd] <- mean(sig_temp[,hd])
  var_T_bias[hd] <- mean(var_T_bias_temp[,hd])
  noise[hd] <- 1/hd_list[hd] * omega2 / (K2*sigma^2) +1
}

plot(mu_hat,sig_hat)

#################### PLOT ####################
data <- data.frame(hd = hd_list, target = (1:length(hd_list)*0)+1, mu_bias = var_mu_bias, sigma_bias = sig_bias, T_bias = var_T_bias, last_noise_bias = noise)

ggplot() +
  geom_line(data=data, aes(x=hd, y=mu_bias, col = "mu_bias"), size = 1) +
  geom_line(data=data, aes(x=hd, y=target, col = "1"), size = 1) +
  geom_line(data=data, aes(x=hd, y=sigma_bias, col = "sigma_bias"), size = 1) +
  geom_line(data=data, aes(x=hd, y=T_bias, col = "T_bias"), size = 1) +
  geom_line(data=data, aes(x=hd, y=last_noise_bias, col = "Last_noise_bias"), size = 1) +
  xlab("Bandwidth, hd") + ylab("Normalized Variance")

print(Sys.time()-p0)
