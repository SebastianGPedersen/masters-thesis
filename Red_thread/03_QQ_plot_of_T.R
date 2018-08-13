setwd(Sys.getenv("masters-thesis"))
#load("../Personal/Fit.Rda") #15-20sek. Loads the fit with Gamma function.
load("Fit.Rdata")
library(ggplot2)
source("Vol/BSS_Sim.R") 
source("Simulation/Noise.R")
source("Estimation/estimates_reloaded.R")

p0 <- Sys.time() #7sec w. 1000 steps, 1000 paths
Nsteps <- 1000
Npaths <- 2

#Create time points
timepoints <- sim.BSS.equidist_times(Nsteps = Nsteps)

#Save covariance matrix
#save.BSS.cov(hVec = timepoints, nPaths = Npaths, S0 = 200, mu_add = 0, type = "Gamma", Fit = Fit)
p0 <- Sys.time() #7sec w. 1000 steps, 1000 paths


#Fit BSS and create lnS
BSS <- sim.BSS(hVec = timepoints, nPaths = Npaths, S0 = 200, mu_add = 0, type = "Gamma", Fit = Fit)
print(Sys.time()-p0)

#save(BSS, file = "Read_thread/BSS.Rda")

#Add noise
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
BSS_w_noise <- sim.addnoise(BSS,omega,rho = 0)


### Parameters for T-calculations
hd <-  5 / (60*24*7*52)
ratio <- 12
lag <- 10

#Time_indices for T-estimator
ten_minutes <- 10 / (6.5*60)
ten_min_index <- ceiling(Nsteps*ten_minutes)
middle_index <- ceiling(Nsteps/2)
ten_last_index <- Nsteps-ten_min_index
desired_indices <- c(ten_min_index, middle_index, ten_last_index)


###Begin T-estimator calculations
#dy
process <- BSS_w_noise
dY <- t(diff(t(process$Y))) #Has to be transposed for 'diff' to work as desired
process$Y <- dY

mu_hat <- est.mu.mat.2.0(process, hd = hd, bandwidth_rescale = T)$mu[,desired_indices]
sig_hat2 <- est.sigma.mat.2.0(process, hv = hd*ratio, lag = lag, bandwidth_rescale = T)$sig[,desired_indices]

T_hat <- sqrt(hd) * mu_hat / sqrt(sig_hat2)



#### Reshape for plot
plot_data_frame <- data.frame(do.call("rbind",
                                      list(cbind(T_hat[,1],rep(" +10min",dim(T_hat)[1])), 
                                           cbind(T_hat[,2],rep(" middle",dim(T_hat)[1])),
                                           cbind(T_hat[,3],rep(" -10min",dim(T_hat)[1])))))

colnames(plot_data_frame) <- c("T_estimator", "Bandwidth")

#Numeric and factor
plot_data_frame$T_estimator <- as.numeric(as.character(plot_data_frame$T_estimator))
plot_data_frame$Bandwidth <- as.factor(plot_data_frame$Bandwidth)


ggplot(plot_data_frame) + 
  stat_qq(aes(sample = T_estimator, colour = Bandwidth)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Quantile of standard normal") + ylab("Quantile of T") +
  theme(legend.position=c(0.87,0.16))

