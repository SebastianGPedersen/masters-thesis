library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 1.6*10^(-5) #What Mathias wrote
omega2 <- omega^2
K2 <- 0.5 #K2

#Constants regarding pre-avg kernel
phi_1 <- 1 #int(g'(x)^2)
phi_2 <- 1/12 #int(g^2)
int_g <- 1/4 #int(g)

#################### PARAMETERS CHANGING WITH N ####################
omega_list <- omega*10^c(-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)

#Initialize list with 5 mean, lower and upper for later plot
all_plot_data <- vector("list", 5)
for (i in 1:5) {
  all_plot_data[[i]]$means <- numeric(length = length(omega_list))
  all_plot_data[[i]]$lower <- all_plot_data[[i]]$means
  all_plot_data[[i]]$upper <- all_plot_data[[i]]$means
}

#################### LOOP OVER N ####################

for (my_omega in 1:length(omega_list)) {
  #my_n <- 1
  #my_n <- length(n_list)
  mat <- 6.5/(24*7*52)#*52*7*24*60*60 #In years
  n <- 23400
  dt <- mat/n #In years
  theta <- 10^(-1)
  k_n <- ceiling(theta*n^(1/2)) #It shouldn't become zero, therefore 'ceiling'
  #hd <- 2*10^3*n^(-1/4) #If miliseconds
  hd <- 10^(-3)*dt^(1/3) #If years
  
  #Check they are reasonable
  n/k_n #We shouldn't use to many in pre-avg
  hd/(k_n*dt) #We should have enough pre-avg within bandwidth
  
  #The index where i calculate T
  desired_index <- floor(n/2) + k_n #burst time + k_n
  
  
  ############ Simulation #########
  Npath <- 300
  settings <- sim.setup(mat=mat, Npath = Npath, Nsteps = n, omega = omega_list[my_omega]) #6.5 hours
  
  Heston <- sim.heston(settings)
  
  # NON NOISE
  Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = 1.5, beta = 0.1, reverse = F)
  Heston_vbdb_noise <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05,c_1 = 0.016, alpha=0.8, reverse = F)
  Heston_jump_noise <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = 0.016, alpha = 0.8)
  
  # OUI NOISE
  Heston_vbdb_large <- list(time = Heston_vbdb_noise$time, Y = Heston_vbdb_noise$X)
  Heston_jump_large <- list(time = Heston_jump_noise$time, Y = Heston_jump_noise$X)
  
  #All paths
  all_paths <- list(Heston, Heston_vbdb_noise, Heston_jump_noise, Heston_vbdb_large, Heston_jump_large)
  
  for (j in 1:length(all_paths)) {
    #j <- 1
    path <- all_paths[[j]]
    
    
    ####### Pre-averaging #########
    #We need to transpose Y for dy to work properly
    Y <- t(as.matrix(path$Y))
    dy <- diff(Y)
    
    # CASES WITH PRE-AVR
    #if(j == 2 || j == 3){
    if(F){
      #Pre-avg (can't take vector)
      pre_y <- matrix(NA, nrow = Npath, ncol = n+1) #includes 0 and n to make same size as Heston$Y
      
      for (i in 1:Npath){
        pre_y[i,] <- c(rep(0,k_n+1),est.NewPreAverage(dy[,i],k_n))
      }
      #Plug back into Heston-data
      path$raw <- path$Y
      path$Y <- pre_y
    }
    else{
      #Plug back into Heston-data
      path$raw <- path$Y
      path$Y <- t(dy)
    }
    
    ######## CALCULATE T estimator ##########
    T_hat <- numeric(length = Npath)
    
    #if(j == 2 || j == 3){
    if(F){
      for (i in 1:Npath){
        #i <- 1
        single_path <- list(Y = path$Y[i,], time = Heston$time, raw = Heston$raw[i,])
        mu_hat <- est.mu.new(data = single_path,hd,t.index = c(desired_index-1,desired_index), kn = k_n)$mu[2]
        sigma_hat_2 <- est.sigma.new(data = single_path,hd, t.index = desired_index, kn = k_n, noisefun = est.noise.iid, theta = theta)$sig
        T_hat[i] <- sqrt(hd/K2) * mu_hat/sqrt(sigma_hat_2) #Correct
      }
    }
    else{
      for (i in 1:Npath){
        #i <- 1
        single_path <- list(Y = path$Y[i,], time = path$time, raw = path$raw[i,])
        mu_hat <- est.mu(data = single_path, hd,t.index = desired_index)$mu[1]
        sigma_hat_2 <- est.sigma.raw(data = single_path, hd,t.index = desired_index)$sig[1]
        T_hat[i] <- sqrt(hd/K2) * mu_hat/sqrt(sigma_hat_2) #Correct
      }
    }
    
    T_hat_clean <- na.omit(T_hat) #Problemer m. negativ sigma nogle gange. Midlertidig løsning.
    
    
    ######## SAVE MEAN AND VARIANCE FOR PLOT #######
    print((length(T_hat)-length(T_hat_clean))/length(T_hat)) #Percentage with na's
    
    all_plot_data[[j]]$means[my_omega] <- mean(T_hat_clean)
    all_plot_data[[j]]$lower[my_omega] <- all_plot_data[[j]]$means[my_omega] -sqrt(var(T_hat_clean))
    all_plot_data[[j]]$upper[my_omega] <- all_plot_data[[j]]$means[my_omega] +sqrt(var(T_hat_clean))
  }
}


#################### PLOT ####################
#Re-shape to data.frame(x, lower, mean, upper, farve)

#Create numbers for colors
all_plot_data[[1]]$factor <- rep("Heston, noise", length(all_plot_data[[1]]$mean))
all_plot_data[[2]]$factor <- rep("Drift burst, noise", length(all_plot_data[[2]]$mean))
all_plot_data[[3]]$factor <- rep("Jump, noise", length(all_plot_data[[3]]$mean))
all_plot_data[[4]]$factor <- rep("Drift burst, large", length(all_plot_data[[4]]$mean))
all_plot_data[[5]]$factor <- rep("Jump, large", length(all_plot_data[[5]]$mean))


#Create a single data_frame
plot_data_frame <- data.frame(omega = omega_list ,
                              lower= all_plot_data[[1]]$lower,
                              mean= all_plot_data[[1]]$mean,
                              upper= all_plot_data[[1]]$upper,
                              factor= all_plot_data[[1]]$factor)

for (i in 2:length(all_plot_data)){
  new_data_frame <- data.frame(omega = omega_list,
                               lower= all_plot_data[[i]]$lower,
                               mean= all_plot_data[[i]]$mean,
                               upper= all_plot_data[[i]]$upper,
                               factor= all_plot_data[[i]]$factor)
  plot_data_frame <- rbind(plot_data_frame,new_data_frame)
}

##### PLOT #####
#qplot(n, mean, data = plot_data_frame, geom = "line", color = factor)

qplot(omega, mean, data = plot_data_frame, geom = "line", color = factor, log = "x") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor), alpha = 0.3)

#if(0 > 1){
#saveRDS(plot_data_frame, file="temp/TdriftVjump.Rda")
#}
#test<-readRDS("temp/TdriftVjump.Rda")
# NOTIFY WHEN COMPLETE
#print(Sys.time())
#shell.exec("https://www.youtube.com/embed/rrVDATvUitA?autoplay=1")
