library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 0.0000225
omega2 <- omega^2
K2 <- 0.5 #K2
theta <- 1/2

#Constants regarding pre-avg kernel
phi_1 <- 1 #int(g'(x)^2)
phi_2 <- 1/12 #int(g^2)
int_g <- 1/4 #int(g)


#################### PARAMETERS CHANGING WITH N ####################
n_list <- c(50, 100, 200, 400, 800, 1600, 2000, 3000, 5000, 7500, 10000, 20000)
means <- numeric(length = length(n_list))
lower <- means
upper <- means



#################### LOOP OVER N ####################

for (my_n in 1:length(n_list)) {
  #my_n <- 1
  #my_n <- length(n_list)
  mat <- 6.5/(24*7*52)*52*7*24*60*60 #In years
  n <- n_list[my_n]
  dt <- mat/n #In years
  k_n <- floor(theta*n^(1/2))
  hd <- 2*10^3*n^(-1/4) #If miliseconds
  #hd <- 10^(-3)*dt^(1/4) #If years
  
  
  ############ Simulation #########
  Npath <- 1000
  settings <- sim.setup(mat=mat, Npath = Npath, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  
  
  ####### Pre-averaging #########
  #We need to transpose Y for dy to work properly
  Y <- t(as.matrix(Heston$Y))
  dy <- diff(Y)
  
  #Pre-avg (can't take vector)
  pre_y <- matrix(NA, nrow = Npath, ncol = n+1) #includes 0 and n to make same size as Heston$Y
  
  for (i in 1:Npath){
    pre_y[i,] <- c(rep(0,k_n+1),est.NewPreAverage(dy[,i],k_n))
  }
    
  #Plug back into Heston-data
  Heston$raw <- Heston$Y
  Heston$Y <- pre_y
  
  
  ######## CALCULATE T estimator ##########
  T_hat <- numeric(length = Npath)
  
  for (i in 1:Npath){
    #i <- 1
    single_path <- list(Y = Heston$Y[i,], time = Heston$time, raw = Heston$raw[i,])
    mu_hat <- est.mu.new(data = single_path,hd,t.index = c(n-1,n), kn = k_n)$mu[2]
    sigma_hat_2 <- est.sigma.new(data = single_path,hd, t.index = n, kn = k_n,noisefun = est.noise.iid, theta = theta)$sig
    T_hat[i] <- sqrt(hd/K2) * mu_hat/sqrt(sigma_hat_2) #Correct
    #T_hat[i] <- sqrt(hd/K2)*mu_hat #mu-estimator. Rimelig konstant varians
    #T_hat[i] <- sigma_hat_2 #sigma_estimator. Stiger ved noise men constant varians uden. Den burde falde
    
  }
  
  T_hat_clean <- na.omit(T_hat) #Problemer m. negativ sigma nogle gange. Midlertidig løsning.
  
  
  ######## SAVE MEAN AND VARIANCE FOR PLOT #######
  print((length(T_hat)-length(T_hat_clean))/length(T_hat)) #Percentage with na's
  
  means[my_n] <- mean(T_hat_clean)
  lower[my_n] <- means[my_n] -sqrt(var(T_hat_clean))
  upper[my_n] <- means[my_n] +sqrt(var(T_hat_clean))
}


#################### PLOT ####################
means
indexes <- 1:length(n_list) #don't take all if na's
plot_data <- data.frame(n = n_list[indexes], means = means[indexes], lower = lower[indexes], upper = upper[indexes])

qplot(n, means, data = plot_data, geom = "line", color = 1)

qplot(n, means, data = plot_data, geom = "line", color = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = 1), alpha = 0.3)

