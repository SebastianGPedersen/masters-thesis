set.seed(100)
library(ggplot2)
library(grid)
library(gridExtra)

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
n <- 23400 /7 #Cheating again
mat <- 6.5/(24*7*52) /7
dt <- mat/n
Npaths <- 1000
sigma2 <- 0.0457 /25
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66

#Because of lack of memory, it is done in loops
n_loops <- 1
desired_index <- n-1

#List to final values
hd <- 5 / (60*24*7*52) #5 min
lag_list <- c(5,10,100)
T_estimator <- matrix(nrow = Npaths,ncol = length(lag_list))

for (memory in 1:n_loops) {
  #memory <- 1
  temp_paths <- Npaths / n_loops
  
  #Heston simulations
  settings <- sim.setup(mat=mat, Npath = temp_paths, Nsteps = n, omega = omega) #6.5 hours
  Heston <- sim.heston(settings)
  
  #### Calculate dY ####
  dY <- t(diff(t(Heston$Y))) #Has to be transposed for 'diff' to work as desired
  data <- list(Y = dY, time = Heston$time)    
  
  for (lag in 1:length(lag_list)) {
    print(paste0("memory = ",memory, ", hd = ",hd,sep = ""))
    sig_hat <- est.sigma.mat(data, hv = hd, t.index = desired_index, lag = lag_list[lag])$sig
    mu_hat <- est.mu.mat(data, hd = hd, t.index = desired_index)$mu
    T_estimator[((memory-1)*temp_paths+1):(memory*temp_paths),lag] <- sqrt(hd) * mu_hat / sqrt(sig_hat)    
  }

}


#################### PLOT ####################
#Re-shape to data-frame
plot_data_frame <- data.frame(do.call("rbind",
                           list(cbind(T_estimator[,1],rep("lag =  5",dim(T_estimator)[1])), 
                                cbind(T_estimator[,2],rep("lag = 10",dim(T_estimator)[1])),
                                cbind(T_estimator[,3],rep("lag = 100",dim(T_estimator)[1])))))

colnames(plot_data_frame) <- c("T_estimator", "Lags")

#Numeric and factor
plot_data_frame$T_estimator <- as.numeric(as.character(plot_data_frame$T_estimator))
plot_data_frame$Lags <- as.factor(plot_data_frame$Lags)


g1 <- ggplot(plot_data_frame, aes(x = T_estimator, fill = Lags)) + 
  xlab("T value") + ylab("Density") + 
  geom_density(adjust = 1.5, alpha = 0.3)+
  theme(legend.position="none")
g1
g2 <- ggplot(plot_data_frame) + 
  stat_qq(aes(sample = T_estimator, colour = Lags)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Quantile of standard normal") + ylab("Quantile of T") +
  theme(legend.position=c(0.87,0.16))
g2

grid.arrange(g1,g2,nrow = 1,
             top = textGrob("Distribution of T-estimator",gp=gpar(fontsize=20)))
