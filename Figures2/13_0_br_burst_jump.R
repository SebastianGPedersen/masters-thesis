library(ggplot2)
library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/add_all.R")
source("Estimation/estimates_reloaded.R")
source("Estimation/estimates_revolution.R")
source("Estimation/rho.R")
source("estimation/teststat.R")

set.seed(100)

####### ESTIMATION PARAMETERS
heston_params <- sim.setup()
lag <- 10


#################### PARAMETERS ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66

#Burn 5min in mu in:
one_min <- 1 / (52*7*24*60)
n_burn <- one_min / dt
t.index <- seq(from = n_burn, to = n, by = 5) #Burn a volatility bandwidth (note 10 in Christensen)

#Burst settings1
c_1_func <- function(alpha) {
  c_1 <- 0 #If alpha == 0
  if (alpha == 0.55) {
    c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)
  } 
  else if (alpha == 0.65){
    c_1 <- (1-alpha)*0.01/(10/(60*24*7*52))^(1-alpha)
  } 
  else if (alpha == 0.75) {
    c_1 <- (1-alpha)*0.015/(10/(60*24*7*52))^(1-alpha)
  }
  return(c_1)
}
c_2_func <- function(beta) {
  c_2 <- 0 #if beta = 0
  if (beta == 0.1) {
    c_2 <- sqrt((1-2*beta)*(8*0.00093*0.25)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.2) {
    c_2 <- sqrt((1-2*beta)*(8*0.00093*0.5)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.3) {
    c_2 <- sqrt((1-2*beta)*(8*0.00093*0.75)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.4) {
    c_2 <- sqrt((1-2*beta)*(8*0.00093*1)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  return(c_2)
}

burstsettings <- list()
burstsettings[[1]] <- sim.burstsetting(jump = F,
                                       alpha = 0,
                                       c_1 = 0,
                                       beta = 0,
                                       c_2 = 0)
burstsettings[[2]] <- sim.burstsetting(jump = F,
                                       alpha = 0.55,
                                       c_1 = c_1_func(0.55),
                                       beta = 0,
                                       c_2 = 0)
burstsettings[[3]] <- sim.burstsetting(jump = T,
                                       alpha = 0.55,
                                       c_1 = c_1_func(0.55),
                                       beta = 0,
                                       c_2 = 0)
burstsettings[[4]] <- sim.burstsetting(jump = F,
                                       alpha = 0,
                                       c_1 = 0,
                                       beta = 0.2,
                                       c_2 = c_2_func(0.2))


#The ratio parameters
ratio_list <- 1:30
h_mu <- 5 / (52*7*24*60)

#Because of lack of memory, it is done in loops
Npaths <- 1000 #15min
n_loops <- ceiling(Npaths/100) #After 50 it just scales linearly if not slower
output_list <- list()
N <- ceiling(Npaths/n_loops)

#Initialize lists used later
Tstar <- numeric(N);  rho<-numeric(N);
temp <- matrix(nrow = length(burstsettings), ncol = length(ratio_list))
output_list <- rep(list(temp),n_loops)

p0 <- Sys.time()

#################### LOOP OVER N ####################
for (memory in 1:n_loops) {
  #memory <- 1
  
  ### Keep track
  print(paste("Loop",memory, "out of ",n_loops))
  print(paste("Expected time left:", 3*round(Npaths-(memory-1)*Npaths/n_loops,0),"seconds")) #One second per path
  
  ### SIMULATE ALL PATHS
  heston_params$Npath <- N
  Heston <- sim.heston(heston_params)
  
  all_simulations <- sim.add_all(Heston = Heston, burst_args = burstsettings)
  
  for (process in 1:length(all_simulations)) {
    #Create dy
    all_simulations[[process]]$Y <- t(diff(t(as.matrix(all_simulations[[process]]$Y))))
  }
  
  for (process in 1:length(all_simulations)) {
    #process <- 1
    mu_hat <- est.mu.mat.2.0(data = all_simulations[[process]], hd = h_mu, bandwidth_rescale = T)$mu[,t.index]
    
    for (ratio in 1:length(ratio_list)) {
      #ratio <- 5
        
      #Estimate T
      sigma2_hat <- est.sigma.mat.3.0(data = all_simulations[[process]], hv = ratio_list[ratio]*h_mu, bandwidth_rescale = T)$sig[,t.index]
      T_hat <- sqrt(h_mu) * mu_hat/ sqrt(sigma2_hat)
      
      #Calculate T_max and thresholds
      for(subpath in 1:N){
        Tstar[subpath] <- max(abs(T_hat[subpath,]))
        # fit rho
        rho[subpath] <- est.rho.MLE(T_hat[subpath,])
      }

      m <- rep(ncol(T_hat),N)
      z<-est.z_quantile(m, rho, 0.95)$qZm
      
      #Save percentage rejection in list
      temp[process,ratio] <- mean(Tstar>=z)
    }
  }
  output_list[[memory]] <- temp
}

print(Sys.time()-p0)

###Take mean over memory loops
output_mean <- matrix(0,nrow = length(all_simulations), ncol = length(ratio_list))

for (memory in 1:n_loops) {
  output_mean <- output_mean + output_list[[memory]]
}
output_mean <- output_mean/n_loops*100

processes <- c("Heston", "Drift burst", "Jump", "Vol burst")

#Create a single data_frame
plot_data_frame <- data.frame(ratios = ratio_list ,
                              value = rep(5,length(ratio_list)),
                              process = "5%")

for (i in 1:nrow(output_mean)){
  new_data_frame <- data.frame(ratios = ratio_list,
                               value = output_mean[i,],
                               process = processes[i])
  plot_data_frame <- rbind(plot_data_frame,new_data_frame)
}


##### PLOT #####
ggplot(plot_data_frame, aes(ratios, value, color = process)) +
  geom_line() +
  xlab(TeX("Bandwidth ratio, C")) + ylab(TeX('$P(T^* > q_{95})$')) + # 'Rejection percentage') +
  ggtitle("Rejection of T-estimator") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
  
#save(plot_data_frame, file="Figures2/Saved_data_for_plots/13_0_br_burst_jump.Rda")

