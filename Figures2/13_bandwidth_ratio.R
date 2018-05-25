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


#################### PARAMETERS ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 80
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
lag <- 10
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66
h_mu <- 5 / (60*24*7*52) #5min bandwidth of drift-estimator
T_interval <- 5 / (60*60*24*7*52) #5sec between T_calculations
threshold <- q95(mat/T_interval)

#The ratio parameters
ratio_list <- seq(1,20,by = 1) #From 1 to 20

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
#n_loops <- floor(Npaths/100) #100 paths in every loop
n_loops <- floor(Npaths/8) #10 loops with 8 in every

#Initialize lists to hold final values
heston <- matrix(nrow = n_loops, ncol = length(ratio_list))
jump <- heston
small_burst <- heston
large_burst <- heston

rejection_list <- list(heston,jump,small_burst,large_burst)


#################### LOOP OVER N ####################
for (memory in 1:n_loops) {
  #memory <- 1
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)

  #The index where i calculate the T's
  desired_indices <- floor(T_interval/dt)*(1:(mat/T_interval-1))
  desired_indices <- desired_indices[desired_indices >  max(ratio_list)*h_mu/dt] #Burn max(sigma_bandwidth) in

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

  
  #All paths
  all_paths <- list(Heston, Heston_jump, Heston_small_vbdb, Heston_large_vbdb) #Should be same order as rejection_list

  #Calculate T for every path
  for (j in 1:length(all_paths)) {
    #j <- 1
    path <- all_paths[[j]]
    #We need to transpose Y for dy to work properly
    Y <- t(as.matrix(path$Y))
    dy <- diff(Y)
    
    #Plug bag into path-list
    path$Y <- t(dy)
    
    ######## CALCULATE T estimator ##########
    for (ratio_index in 1:length(ratio_list)) {
      #ratio_index <- 1
      p0 <- Sys.time()
      print(paste0("memory = ",memory, ", path = ",j, ", ratio_index = ", ratio_index, sep = ""))
      
      mu_hat <- est.mu.mat.next(data = path, h_mu,t.index = desired_indices)$mu
      sigma_hat2 <- est.sigma.mat.next(data = path, h_mu*ratio_list[ratio_index],t.index = desired_indices)$sig
      T_hat <- max(abs(sqrt(h_mu)*mu_hat/sqrt(sigma_hat2)))
      
      rejection_list[[j]][memory,ratio_index] <- sum((T_hat > threshold)) / length(T_hat) #Save rejection percentage
      print((Sys.time()-p0)*length(ratio_list)*length(all_paths)*n_loops) #20 timer med kun lag = 10 og Npaths = 80
    }
  }
}


#################### PLOT ####################
#Re-shape to data.frame(x, lower, mean, upper, farve)

all_paths <- list(Heston, Heston_vb, Heston_vbdb, Heston_jump, Heston_vbjump)

#Create numbers for colors
all_plot_data[[1]]$process <- rep("Heston", length(all_plot_data[[1]]$mean))
all_plot_data[[2]]$process <- rep("+ volatility burst", length(all_plot_data[[2]]$mean))
all_plot_data[[3]]$process <- rep("+ drift burst & volatility burst", length(all_plot_data[[3]]$mean))
all_plot_data[[4]]$process <- rep("+ jump", length(all_plot_data[[4]]$mean))
all_plot_data[[5]]$process <- rep("+ jump & volatility burst", length(all_plot_data[[5]]$mean))


#Create a single data_frame
plot_data_frame <- data.frame(n = n_list ,
                              lower= all_plot_data[[1]]$lower,
                              mean= all_plot_data[[1]]$mean,
                              upper= all_plot_data[[1]]$upper,
                              process= all_plot_data[[1]]$process)

for (i in 2:length(all_plot_data)){
  new_data_frame <- data.frame(n = n_list,
                               lower= all_plot_data[[i]]$lower,
                               mean= all_plot_data[[i]]$mean,
                               upper= all_plot_data[[i]]$upper,
                               process= all_plot_data[[i]]$process)
  plot_data_frame <- rbind(plot_data_frame,new_data_frame)
}

##### PLOT #####
qplot(n, mean, data = plot_data_frame, geom = "line", color = process) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = process), alpha = 0.3) +
  xlab("Number of observations") + ylab(TeX('$ T-estimator \\pm sd$'))

#Save dataframe for later
saveRDS(plot_data_frame, file="Figures2/Saved_data_for_plots/04_T-estimator2_without_noise.Rda")

print(Sys.time()-p0) #approx 10 min with max(n) = 60k and npaths = 500
