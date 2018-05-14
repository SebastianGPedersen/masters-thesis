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
Npaths <- 400
sigma <- 0.0225
bandwidth_ratio <- 6 #For a normal distribution
(noise_ratio <- omega/sigma*sqrt(n))

#Alternative
#noise_ratio <- 1
#omega <- noise_ratio*sigma/sqrt(n)

#Lag length choice with eq. (3) in Barndorff-Nielsen (2009)
c_star <- 3.5134
xi <- sqrt(omega^2/sqrt(mat^2*sigma^4))
lag <- c_star*xi^(4/5)*n^(3/5)

#Because of lack of memory, it is done in loops
n_loops <- 4

#List to final values
h_list <- c(2,5,10) / (60*24*7*52) #2min, 5min and 10min
T_estimator <- matrix(nrow = Npaths,ncol = length(h_list))


for (memory in 1:n_loops) {
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  set.seed(1000*memory)
  
  #Heston simulations
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  
  for (h_mu in 1:length(h_list)) {
    #h_mu <- 1
    desired_index <- n-1 #Just takes last index so K has as many obs as possible
    
    for (i in 1:temp_paths) {
      #i <- 1
      #Heston
      single_path <- list(Y = diff(Heston$Y[i,]), time = Heston$time)
      sig_hat <- est.sigma(single_path, hv = h_list[h_mu]*bandwidth_ratio, t.index = desired_index, lag = lag)$sig[1]
      mu_hat <- est.mu(data = single_path, hd = h_list[h_mu], t.index = desired_index)$mu[1]
      T_estimator[(memory-1)*temp_paths+i,h_mu] <- sqrt(h_list[h_mu])*mu_hat/sqrt(sig_hat) #K2*sigma cancels out and we end with sqrt(h_n)*\hat{\mu} / sqrt(\hat{\Sigma}^2)
    }
  }
}


#################### PLOT ####################
#Re-shape to data-frame
plot_data_frame <- data.frame(do.call("rbind",
                           list(cbind(T_estimator[,1],rep(" 2min",dim(T_estimator)[1])), 
                                cbind(T_estimator[,2],rep(" 5min",dim(T_estimator)[1])),
                                cbind(T_estimator[,3],rep("10min",dim(T_estimator)[1])))))

colnames(plot_data_frame) <- c("T_estimator", "Bandwidth")

#Numeric and factor
plot_data_frame$T_estimator <- as.numeric(as.character(plot_data_frame$T_estimator))
plot_data_frame$Bandwidth <- as.factor(plot_data_frame$Bandwidth)

ggplot(plot_data_frame) + 
  stat_qq(aes(sample = T_estimator, colour = Bandwidth)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Quantile of standard normal") + ylab("Quantile of T")
