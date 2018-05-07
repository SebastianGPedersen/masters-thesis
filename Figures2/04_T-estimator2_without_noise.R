library(ggplot2)
library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

p0 <- Sys.time()

#seed
set.seed(100)

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 0
K2 <- 0.5 #K2

#Burst settings
alpha <- 0.8
beta <- 0.1
c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)
c_2 <- sqrt((1-2*beta)*0.001^2/(10/(60*24*7*52))^(1-2*beta))

#################### PARAMETERS CHANGING WITH N ####################
#n = 60k and npaths = 500 is absolute max that my computer can keep in memory
n_list <- c(50, 100, 200, 400, 800, 1600, 2000, 3000, 5000, 7500, 10000, 20000, 30000, 40000, 60000)

#Initialize list with 5 mean, lower and upper for later plot
n_processes <- 5

all_plot_data <- vector("list", n_processes)
for (i in 1:n_processes) {
  all_plot_data[[i]]$means <- numeric(length = length(n_list))
  all_plot_data[[i]]$lower <- all_plot_data[[i]]$means
  all_plot_data[[i]]$upper <- all_plot_data[[i]]$means
}

#################### LOOP OVER N ####################

for (my_n in 1:length(n_list)) {
  #my_n <- 1
  print(my_n)
  #my_n <- length(n_list)
  mat <- 6.5/(24*7*52)#*52*7*24*60*60 #In years
  n <- n_list[my_n]
  dt <- mat/n #In years
  hd <- 10^(-2)*dt^(1/2) #If years
  
  #Check that they are reasonable
  hd/dt
  
  #The index where i calculate T
  desired_index <- ceiling((n+1)/2) #because it includes t_0 and T, we take (n+1)/2.
  
  
  ############ Simulation #########
  Npath <- 500
  settings <- sim.setup(mat=mat, Npath = Npath, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = c_2, beta = beta)
  Heston_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05, c_1 = c_1, alpha = alpha)
  Heston_jump <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = c_1, alpha = alpha)
  Heston_vbjump <- sim.addjump(Heston_vb, burst_time = 0.5, interval_length = 0.05, c_1 = c_1, alpha = alpha)
  
  #All paths
  all_paths <- list(Heston, Heston_vb, Heston_vbdb, Heston_jump, Heston_vbjump)

  for (j in 1:length(all_paths)) {
    #j <- 1
    path <- all_paths[[j]]
    
    #We need to transpose Y for dy to work properly
    Y <- t(as.matrix(path$Y))
    dy <- diff(Y)
    
    #Plug bag into path-list
    path$Y <- t(dy)
    
    ######## CALCULATE T estimator ##########
    T_hat <- numeric(length = Npath)
    
    for (i in 1:Npath){
      #i <- 1
      single_path <- list(Y = path$Y[i,], time = path$time)
      mu_hat <- est.mu(data = single_path, hd,t.index = desired_index)$mu[1]
      sigma_hat_2 <- est.sigma.raw(data = single_path, hd,t.index = desired_index)$sig[1]
      T_hat[i] <- sqrt(hd/K2) * mu_hat/sqrt(sigma_hat_2)
    }
    
    T_hat_clean <- na.omit(T_hat) #There MIGHT be negtive sigma-hat but it is highly unlikely
    
    ######## SAVE MEAN AND VARIANCE FOR PLOT #######

    all_plot_data[[j]]$means[my_n] <- mean(T_hat_clean)
    all_plot_data[[j]]$lower[my_n] <- all_plot_data[[j]]$means[my_n] -sqrt(var(T_hat_clean))
    all_plot_data[[j]]$upper[my_n] <- all_plot_data[[j]]$means[my_n] +sqrt(var(T_hat_clean))
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
saveRDS(plot_data_frame, file="Figures2/Saved_data_for_plots/04_T-estimator2_wihtout_noise.Rda")

print(Sys.time()-p0) #approx 10 min with max(n) = 60k and npaths = 500
