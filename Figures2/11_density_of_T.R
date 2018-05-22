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
omega2 <- 2.64*10^(-10)*25
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400 /7 #I cheat by dividing n and mat with 7 (no influence, because same dt and I take last index, so the first 6/7 obs have no influence anyway)
mat <- 6.5/(24*7*52) /7
dt <- mat/n
Npaths <- 1000
sigma2 <- 0.0457/25
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sigma*sqrt(dt)))
lag <- 100

#Because of lack of memory, it is done in loops
n_loops <- 10

#List to final values
hd <- 1 / (60*24*7*52)
T_estimator <- numeric(length = Npaths)

for (memory in 1:n_loops) {
  #memory <- 1
  print(paste("memory = ",memory, sep = ""))
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  
  #Heston simulations
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)

  #### Calculate dY ####
  dY <- t(diff(t(Heston$Y))) #Has to be transposed for 'diff' to work as desired
  data <- list(Y = dY, time = Heston$time)    
  
  desired_index <- n-1
  sig_hat <- est.sigma.mat(data, hv = hd, t.index = desired_index, lag = lag)$sig
  mu_hat <- est.mu.mat(data, hd = hd, t.index = desired_index)$mu
  T_estimator[((memory-1)*temp_paths+1):(memory*temp_paths)] <- sqrt(hd) * mu_hat / sqrt(sig_hat)

}


#################### PLOT ####################
#Re-shape to data-frame
plot_data_frame <- data.frame(T_estimator)
ggplot(data = plot_data_frame) + 
  aes_string(T_estimator) + 
  xlab("T value") + ylab("Density") + 
  geom_density(adjust = 1, fill = gray(0.5)) +
  ggtitle("Density of T-estimator") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
