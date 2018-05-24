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
lag <- 10

#Because of lack of memory, it is done in 10 loops
n_loops <- 5

#List to final values
h_mu <- 5 / (60*24*7*52)
ratio_list <- seq(2,20,by = 1) / 2

#Initialize lists to hold final values
heston <- matrix(nrow = n_loops, ncol = length(ratio_list))
jump <- matrix(nrow = n_loops, ncol = length(ratio_list))

rejection_list <- list(heston, jump)


for (memory in 1:n_loops) {
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  
  #Heston simulations
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  Heston_jump_small <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = 0.3, alpha = 0.55)

  all_paths <- list(Heston, Heston_jump_small) #Should be same order as rejection_list
  
  #################### LOOP OVER PATHS AND H_D ####################
  
  for (j in 1:length(all_paths)) {
    #j <- 1
    path <- all_paths[[j]]
    
    #Transform Y to dY in path$Y
    path$Y <- t(diff(t(as.matrix(path$Y))))
    
    #Extract the relevant part
    path$Y <- path$Y[,1:desired_index]
    path$time <- path$time[1:desired_index]
    path$X <- NULL
    path$vol <- NULL
    
    #H_d from 1min to 10min
  
    for (hd in 1:length(ratio_list)) {
      
      #hd <- 1
      desired_index <- n/2 #
  
      print(paste0("memory = ",memory, ", path = ",j, ", ratio_index = ", hd, sep = ""))
      
      mu_hat <- est.mu.mat(data = path, hd = h_mu,t.index = desired_index)$mu
      sigma_hat2 <- est.sigma.mat(data = path, hv = h_mu*ratio_list[hd],t.index = desired_index,lag = lag)$sig
      T_hat <- abs(sqrt(h_mu)*mu_hat/sqrt(sigma_hat2))
      
      rejection_list[[j]][memory,hd] <- sum((T_hat > threshold)) / length(T_hat) #Save rejection percentage
      #print((Sys.time()-p0)*length(ratio_list)*length(all_paths)*n_loops)
    }
  }
}

#Save mean(T_above_percentage)
Heston_above <- numeric(length = length(ratio_list))
jump_above <- Heston_above

for (hd in 1:length(ratio_list)) {
  Heston_above[hd] <- mean(rejection_list[[1]][,hd])
  jump_above[hd] <- mean(rejection_list[[2]][,hd])
}


#################### PLOT THE PERCENTAGE ABOVE ####################
data <- data.frame(hd = ratio_list, 
                   target = (1:length(ratio_list)*0)+0.05, 
                   Heston_above = Heston_above,
                   jump_above = jump_above)

ggplot() +
  geom_line(data=data, aes(x=hd, y=target, col = "5 percent"), size = 1) +
  geom_line(data=data, aes(x=hd, y=Heston_above, col = "Heston"), size = 1) +
  geom_line(data=data, aes(x=hd, y=jump_above, col = "Jump"), size = 1) +
  xlab("Ratio of h_sigma / h_mu") + ylab("Percentage above T_95") +
  ggtitle("Rejection rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

print(Sys.time()-p0)

