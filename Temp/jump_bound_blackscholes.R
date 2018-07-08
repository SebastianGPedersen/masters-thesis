library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/BlackScholes.R")
source("Simulation/Jumps.R")
source("Simulation/Heston.R")
source("Estimation/estimates_reloaded.R")

#seed
set.seed(100)

i <- desired_indices
#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 0
K2 <- 0.5 #K2
nsteps <- 23400
desired_indices <- ceiling((nsteps-1)/2)
threshold <- qnorm(0.975)
h_mu <- 5/(52*7*24*60)
  
### Memory loop
Npaths <- 500
#n_loops <- ceiling(Npaths/100)
n_loops <- 1
  
#Heston settings
settings <- sim.setup(mat=6.5/(24*7*52), Npath = ceiling(Npaths/n_loops), Nsteps = nsteps, omega = 0) #6.5 hours

#Burst settings
q <- threshold
c <- sqrt(2)
sigma <- sqrt(settings$theta)
k_max <- c*sigma/sqrt(q^2-c^2)
J <- k_max*sqrt(h_mu) #0.06%

alpha <- 0.8
c_1 <- (1-alpha)*J/(10/(60*24*7*52))^(1-alpha)

sd(BS$Y[1,1:(ncol(BS$Y)-1)]-BS$Y[1,2:ncol(BS$Y)])

#################### LOOP OVER N ####################
rejection <- numeric(length = n_loops)
for (memory in 1:n_loops) {
  #memory <- 1
  print(paste("Loop",memory, "out of ",n_loops))
  
  #Heston
  BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = 0, Nsteps = nsteps, Npath = ceiling(Npaths/n_loops))
  BS_jump <- sim.addjump(BS, burst_time = 0.5, interval_length = 0.05, c_1 = c_1, alpha = alpha)
  path <- BS_jump
  #path <- BS
  
  #Create dy
  Y <- t(as.matrix(path$Y))
  dy <- diff(Y)
  path$Y <- t(dy)
  
  #T-estimator
  mu_hat <- est.mu.mat.2.0(data = path, h_mu)$mu[,desired_indices-1]
  sigma_hat_2 <- est.sigma.raw.mat.2.0(data = path, h_mu)$sig[,desired_indices-1]
  T_hat <- sqrt(h_mu/K2) * mu_hat/sqrt(sigma_hat_2)
  
  T_hat[is.na(T_hat)] <- 0 #There MIGHT be negtive sigma-hat but it is highly unlikely
  
  rejection[memory] <- mean((T_hat < -threshold))
}
rejection_perc <- mean(rejection)*100

### Remember jump is in fact smaller than J
### Check distribution results of mu
mu_dist <- sqrt(h_mu)*mu_hat

#Check mean. Mean is only half of what it's supposed to
mean(mu_dist)
1/sqrt(h_mu)*(-J)

#Check variance. Seems correct
sd(mu_dist)
sqrt(K2)*sigma

### Check distribution results of sigma
sigma^2+1/h_mu*J^2
mean(sigma_hat_2)
sd(sigma_hat_2)


x_min <- -(60*10):(60*30)/60
normal_distribution <- 2.5
max_rejection <- (1-pnorm(sqrt(q^2-c^2)))*100

df <- data.frame('rejection_percentage' = rejection_perc, 'min_after_jump' = x_min,
                 'threshold' = rep(normal_distribution,length(x_min)),
                 'max_rejection' = rep(max_rejection,length(x_min)))
ggplot(df,aes(min_after_jump,rejection_percentage)) +
  geom_smooth()

#geom_point(min_after_jump,rejection_percentage) +
geom_line(min_after_jump,threshold) +
geom_line(min_after_jump,max_rejection)



#################### PLOT ####################
#Re-shape to data.frame(x, lower, mean, upper, farve)

all_paths <- list(Heston, Heston_vb, Heston_vbdb, Heston_jump, Heston_vbjump)

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
qplot(n, mean, data = plot_data_frame, geom = "line", color = process) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = process), alpha = 0.3) +
  xlab("Number of observations") + ylab(TeX('$ T-estimator \\pm sd$'))

#Save dataframe for later
save(plot_data_frame, file="Figures2/Saved_data_for_plots/04_T-estimator2_without_noise.Rda")

print(Sys.time()-p0) #approx 10 min with max(n) = 60k and npaths = 500
