library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

#################### PARAMETERS THAT DON'T CHANGE ####################
#omega <- 0.0000225
omega <- 0
#omega2 <- omega^2
K2 <- 0.5 #K2
#theta <- 1/2

#Constants regarding pre-avg kernel
#phi_1 <- 1 #int(g'(x)^2)
#phi_2 <- 1/12 #int(g^2)
#int_g <- 1/4 #int(g)


#################### PARAMETERS CHANGING WITH N ####################
n_list <- c(50, 100, 200, 400, 800, 1600, 2000, 3000, 5000, 7500, 10000, 20000, 30000, 40000, 50000, 60000)

#Initialize list with 5 mean, lower and upper for later plot
all_plot_data <- vector("list", 5)
for (i in 1:5) {
  all_plot_data[[i]]$means <- numeric(length = length(n_list))
  all_plot_data[[i]]$lower <- all_plot_data[[i]]$means
  all_plot_data[[i]]$upper <- all_plot_data[[i]]$means
}

#################### LOOP OVER N ####################

for (my_n in 1:length(n_list)) {
  #my_n <- 1
  #my_n <- length(n_list)
  mat <- 6.5/(24*7*52)#*52*7*24*60*60 #In years
  n <- n_list[my_n]
  dt <- mat/n #In years
  #k_n <- floor(theta*n^(1/2))
  #hd <- 10^(-3)*dt^(1/4) #If years
  hd <- 10^(-1.5)*dt^(1/2) #If years
  
  #Check that they are reasonable
  hd/dt
  #mat/k_n
  
  #The index where i calculate T
  desired_index <- ceiling((n+1)/2)# + k_n #burst time + k_n
  
  
  ############ Simulation #########
  Npath <- 250
  settings <- sim.setup(mat=mat, Npath = Npath, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  
  #alpha < beta + 1/2. Burst + Jump
  Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = 0.03, beta = 0.45)
  Heston_vbdb_small <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05,c_1 = 0.3,alpha=0.55)
  Heston_jump_small <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = 0.3, alpha = 0.55)
  
  #alpha > beta + 1/2. Burst + Jump
  Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = 1.5, beta = 0.1)
  Heston_vbdb_large <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05,c_1 = 0.016,alpha=0.8)
  Heston_jump_large <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = 0.016, alpha = 0.8)
  
  #All paths
  all_paths <- list(Heston, Heston_vbdb_small, Heston_jump_small, Heston_vbdb_large, Heston_jump_large)

  for (j in 1:length(all_paths)) {
    #j <- 1
    path <- all_paths[[j]]
    
    ####### Pre-averaging #########
    #We need to transpose Y for dy to work properly
    Y <- t(as.matrix(path$Y))
    dy <- diff(Y)
    
    #Pre-avg (can't take vector)
    #pre_y <- matrix(NA, nrow = Npath, ncol = n+1) #includes 0 and n to make same size as Heston$Y
    
    #for (i in 1:Npath){
    #  pre_y[i,] <- c(rep(0,k_n+1),est.NewPreAverage(dy[,i],k_n))
    #}
      
    #Plug back into Heston-data
    path$raw <- path$Y
    #path$Y <- pre_y
    
    #Temp
    path$Y <- t(dy)
    
    ######## CALCULATE T estimator ##########
    T_hat <- numeric(length = Npath)
    
    for (i in 1:Npath){
      #i <- 1
      single_path <- list(Y = path$Y[i,], time = path$time, raw = path$raw[i,])
      mu_hat <- est.mu(data = single_path, hd,t.index = desired_index)$mu[1]
      sigma_hat_2 <- est.sigma.raw(data = single_path, hd,t.index = desired_index)$sig[1]
      T_hat[i] <- sqrt(hd/K2) * mu_hat/sqrt(sigma_hat_2) #Correct
      #T_hat[i] <- sqrt(hd/K2)*mu_hat #mu-estimator. Rimelig konstant varians
      #T_hat[i] <- sigma_hat_2 #sigma_estimator. Stiger ved noise men constant varians uden. Den burde falde
    }
    
    T_hat_clean <- na.omit(T_hat) #Problemer m. negativ sigma nogle gange. Midlertidig løsning.
    
    
    ######## SAVE MEAN AND VARIANCE FOR PLOT #######
    print((length(T_hat)-length(T_hat_clean))/length(T_hat)) #Percentage with na's
    
    all_plot_data[[j]]$means[my_n] <- mean(T_hat_clean)
    all_plot_data[[j]]$lower[my_n] <- all_plot_data[[j]]$means[my_n] -sqrt(var(T_hat_clean))
    all_plot_data[[j]]$upper[my_n] <- all_plot_data[[j]]$means[my_n] +sqrt(var(T_hat_clean))
  }
}


#################### PLOT ####################
#Re-shape to data.frame(x, lower, mean, upper, farve)

#Create numbers for colors
all_plot_data[[1]]$factor <- rep("Heston", length(all_plot_data[[1]]$mean))
all_plot_data[[2]]$factor <- rep("Drift burst, small", length(all_plot_data[[2]]$mean))
all_plot_data[[3]]$factor <- rep("Jump, small", length(all_plot_data[[3]]$mean))
all_plot_data[[4]]$factor <- rep("Drift burst, large", length(all_plot_data[[4]]$mean))
all_plot_data[[5]]$factor <- rep("Jump, large", length(all_plot_data[[5]]$mean))


#Create a single data_frame
plot_data_frame <- data.frame(n = n_list ,
                              lower= all_plot_data[[1]]$lower,
                              mean= all_plot_data[[1]]$mean,
                              upper= all_plot_data[[1]]$upper,
                              factor= all_plot_data[[1]]$factor)

for (i in 2:length(all_plot_data)){
  new_data_frame <- data.frame(n = n_list,
                               lower= all_plot_data[[i]]$lower,
                               mean= all_plot_data[[i]]$mean,
                               upper= all_plot_data[[i]]$upper,
                               factor= all_plot_data[[i]]$factor)
  plot_data_frame <- rbind(plot_data_frame,new_data_frame)
}

##### PLOT #####
#qplot(n, mean, data = plot_data_frame, geom = "line", color = factor)

qplot(n, mean, data = plot_data_frame, geom = "line", color = factor) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor), alpha = 0.3)

