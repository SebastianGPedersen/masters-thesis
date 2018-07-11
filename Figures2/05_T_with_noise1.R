library(ggplot2)
library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

#seed
set.seed(100)

p0 <- Sys.time()

#################### PARAMETERS THAT DON'T CHANGE ####################
omega2 <- 2.64*10^(-10) #What Mathias wrote
omega <- sqrt(omega2)
K2 <- 0.5 #K2

#Burst settings
alpha <- 0.55
beta <- 0.45
c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)
c_2 <- sqrt((1-2*beta)*0.001^2/(10/(60*24*7*52))^(1-2*beta))


#Constants regarding pre-avg kernel
phi_1 <- 1 #int(g'(x)^2)
phi_2 <- 1/12 #int(g^2)
int_g <- 1/4 #int(g)


#################### PARAMETERS CHANGING WITH N ####################
n_list <- c(800, 1600, 2000, 3000, 5000, 7500, 10000,
            20000, 30000, 40000, 60000)

#Initialize list with 5 mean, lower and upper for later plot

#Initialize list with 5 mean, lower and upper for later plot
n_processes <- 4

all_plot_data <- vector("list", n_processes)
for (i in 1:n_processes) {
  all_plot_data[[i]]$means <- numeric(length = length(n_list))
  all_plot_data[[i]]$lower <- all_plot_data[[i]]$means
  all_plot_data[[i]]$upper <- all_plot_data[[i]]$means
}

#################### LOOP OVER N ####################

for (my_n in 1:length(n_list)) {
  print(my_n)
  #my_n <- 8
  #my_n <- length(n_list)
  mat <- 6.5/(24*7*52)#*52*7*24*60*60 #In years
  n <- n_list[my_n]
  dt <- mat/n #In years
  theta_1 <- 10^(-1)
  k_n <- ceiling(theta_1*n^(1/2)) #It shouldn't become zero, therefore 'ceiling'
  #hd <- 2*10^3*n^(-1/4) #If miliseconds
  theta_2 <- 10^(-1)
  hd <- mat*theta_2*n^(-1/3) #If years
  
  #Check they are reasonable
  #n/k_n #We shouldn't use to many in pre-avg
  #hd/(k_n*dt) #We should have enough pre-avg within bandwidth
  
  #The index where i calculate T
  desired_index <- ceiling((n+1)/2) + k_n-1 #burst time + k_n - 1
  
  
  ############ Simulation #########
  Npath <- 400
  settings <- sim.setup(mat=mat, Npath = Npath, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = c_2, beta = beta, reverse = F, recenter = F)
  Heston_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05, c_1 = c_1, alpha = alpha, reverse = F)
  Heston_jump <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = c_1, alpha = alpha)

  #All paths
  all_paths <- list(Heston, Heston_vb, Heston_vbdb, Heston_jump)

  for (j in 1:length(all_paths)) {
    #j <- 1
    path <- all_paths[[j]]
    
    ####### Pre-averaging #########
    #We need to transpose Y for dy to work properly
    Y <- t(as.matrix(path$Y))
    dy <- diff(Y)
    
    #Pre-avg (can't take vector)
    pre_y <- matrix(NA, nrow = Npath, ncol = n+1) #includes 0 and n to make same size as Heston$Y
    
    for (i in 1:Npath){
      pre_y[i,] <- c(rep(0,k_n+1),est.NewPreAverage(dy[,i],k_n))
    }
      
    #Plug back into Heston-data
    path$raw <- path$Y
    path$Y <- pre_y
    
    
    ######## CALCULATE T estimator ##########
    T_hat <- numeric(length = Npath)
    
    for (i in 1:Npath){
      #i <- 1
      single_path <- list(Y = path$Y[i,], time = Heston$time, raw = Heston$raw[i,])
      mu_hat <- est.mu.new(data = single_path,hd,t.index = c(desired_index-1,desired_index), kn = k_n)$mu[2]
      sigma_hat_2 <- est.sigma.new(data = single_path,hd, t.index = desired_index, kn = k_n,noisefun = est.noise.iid, theta = theta_1)$sig
      T_hat[i] <- sqrt(hd/K2) * mu_hat/sqrt(sigma_hat_2) #Correct
    }
    
    T_hat_clean <- na.omit(T_hat) #There MIGHT be negtive sigma-hat but it is highly unlikely
    
    ######## SAVE MEAN AND VARIANCE FOR PLOT #######
    all_plot_data[[j]]$means[my_n] <- mean(T_hat_clean)
    all_plot_data[[j]]$lower[my_n] <- all_plot_data[[j]]$means[my_n] - sqrt(var(T_hat_clean))
    all_plot_data[[j]]$upper[my_n] <- all_plot_data[[j]]$means[my_n] +sqrt(var(T_hat_clean))
  }
}


#################### PLOT ####################
#Re-shape to data.frame(x, lower, mean, upper, farve)

#Create numbers for colors
all_plot_data[[1]]$process <- rep("Heston", length(all_plot_data[[1]]$mean))
all_plot_data[[2]]$process <- rep("+ volatility burst", length(all_plot_data[[2]]$mean))
all_plot_data[[3]]$process <- rep("+ drift burst & volatility burst", length(all_plot_data[[3]]$mean))
all_plot_data[[4]]$process <- rep("+ jump", length(all_plot_data[[4]]$mean))


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
ggplot(plot_data_frame, aes(n, mean, color = process)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = process), alpha = 0.3) +
  xlab("Number of observations") + ylab(TeX('$ T-estimator \\pm sd$')) +
  ggtitle("T-estimator for noisy processes") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

#Save dataframe for later
save(plot_data_frame, file="Figures2/Saved_data_for_plots/05_T-estimator1_with_noise.Rda")

print(Sys.time()-p0) #approx 10 min with max(n) = 60k and npaths = 500
