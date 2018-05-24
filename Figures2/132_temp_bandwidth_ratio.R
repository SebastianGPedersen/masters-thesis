library(ggplot2)
library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")
source("Estimation/gumbel.R")
p0 <- Sys.time()

#I cheat with 5 things in this temporary study:
#1) Lag is only 10
#2) I don't evaluate max(T) in every point but only T at jump and drift time.
#3) I evaluate it against Gumbel and not AR(1) (simply because this hasn't been coded yet)
#4) I don't rescale the beginning (simply because this hasn't been coded yet)
#5) I only have 100 paths instead of 1000

#With all the cheat above it only takes a minute
#Without cheating it takes weeks.
#I don't think any of the above should change the results significantly.


#################### PARAMETERS ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 200 #Temporary. Should be 1000
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
lag <- 10 #Temporary. Should be 100
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66
h_mu <- 5 / (60*24*7*52) #5min bandwidth of drift-estimator
#T_interval <- 5 / (60*60*24*7*52) #5sec between T_calculations
T_interval <- 1 / (60*24*7*52) #1min between T_calculations
threshold <- q95(mat/T_interval)
#threshold <- qnorm(0.975)

#The ratio parameters
ratio_list <- seq(1,41,by = 2)

#Burst settings1
alpha_large <- 0.8
beta_large <- 0.1

#Burst settings2
alpha_small <- 0.55
beta_small <- 0.45

#Constants for bursts
c_1_func <- function(alpha) {(1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)}
c_2_func <- function(beta) {sqrt((1-2*beta)*0.001^2/(10/(60*24*7*52))^(1-2*beta))}

#Because of lack of memory, it is done in loops
n_loops <- floor(Npaths/100) #100 paths in every loop
#n_loops <- floor(Npaths/8) #10 loops with 8 in every

#Initialize lists to hold final values
jump <- matrix(nrow = n_loops, ncol = length(ratio_list))
small_burst <- jump
large_burst <- jump

rejection_list <- list(jump,small_burst,large_burst)


#################### LOOP OVER N ####################
for (memory in 1:n_loops) {
  #memory <- 1
  temp_paths <- Npaths / n_loops
  set.seed(100*memory)

  #The index where i calculate the T's
  desired_index <- n/2
  
  ### Simulations
  settings <- sim.setup(mat=mat, Npath = Npaths, Nsteps = n, omega = omega) #6.5 hours
  
  #Heston
  Heston <- sim.heston(settings)
  #Large bursts
  Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = c_2_func(beta_large), beta = beta_large)
  Heston_large_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05, c_1 = c_1_func(alpha_large), alpha = alpha_large)
  #Small burst
  Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = c_2_func(beta_large), beta = beta_large)
  Heston_small_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05, c_1 = c_1_func(alpha_large), alpha = alpha_large)
  
  #Jump
  Heston_jump <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = c_1_func(alpha_small), alpha = alpha_small) #Same jump independent of alpha_small or alpha_large

  
  #List of all paths
  all_paths <- list(Heston_jump, Heston_small_vbdb, Heston_large_vbdb) #Should be same order as rejection_list

  #Calculate T for every path
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
    
    ######## CALCULATE T estimator ##########
    for (ratio_index in 1:length(ratio_list)) {
      #ratio_index <- 1
      #ratio_index <- length(ratio_list)
      p0 <- Sys.time()
      print(paste0("memory = ",memory, ", path = ",j, ", ratio_index = ", ratio_index, sep = ""))
      
      mu_hat <- est.mu.mat(data = path, hd = h_mu,t.index = desired_index)$mu
      sigma_hat2 <- est.sigma.mat(data = path, hv = h_mu*ratio_list[ratio_index],t.index = desired_index,lag = lag)$sig
      T_hat <- abs(sqrt(h_mu)*mu_hat/sqrt(sigma_hat2))
      
      rejection_list[[j]][memory,ratio_index] <- sum((T_hat > threshold)) / length(T_hat) #Save rejection percentage
      print((Sys.time()-p0)*length(ratio_list)*length(all_paths)*n_loops)
    }
  }
}

#Save mean(T_above_percentage)
rejection_mean <- list()
for (i in 1:length(rejection_list)) {rejection_mean[[i]] <- numeric(length = length(ratio_list))}

for (hd in 1:length(ratio_list)) {
  for (x in 1:length(rejection_list)){
    rejection_mean[[x]][hd] <- mean(rejection_list[[x]][,hd])
  }
}

#################### PLOT THE PERCENTAGE ABOVE ####################
data <- data.frame(hd = ratio_list, 
                   target = (1:length(ratio_list)*0)+0.05, 
                   jump = rejection_mean[[1]],
                   small_burst = rejection_mean[[2]],
                   large_burst = rejection_mean[[3]])

ggplot() +
  geom_line(data=data, aes(x=hd, y=target, col = "5 percent"), size = 1) +
  geom_line(data=data, aes(x=hd, y=jump, col = "Jump"), size = 1) +
  geom_line(data=data, aes(x=hd, y=small_burst, col = "Burst"), size = 1) +
  #geom_line(data=data, aes(x=hd, y=large_burst, col = "Large burst"), size = 1) +
  xlab("Ratio of h_sigma / h_mu") + ylab("Percentage above T_95") +
  ggtitle("Rejection rate") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

### Save data

print(Sys.time()-p0)
