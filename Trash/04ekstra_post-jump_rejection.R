library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Jumps.R")
source("Estimation/estimates_reloaded.R")

p0 <- Sys.time()

#seed
set.seed(100)

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 0
K2 <- 0.5 #K2
nsteps <- 23000
desired_indices <- (nsteps-1)/2 + -(60*10):(60*30) #10min before, 30min after
threshold <- qnorm(0.975)
  
#Burst settings
alpha <- 0.8
c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)


### Memory loop
Npaths <- 1000
n_loops <- ceiling(Npaths/100)

#Heston settings
settings <- sim.setup(mat=6.5/(24*7*52), Npath = ceiling(Npaths/n_loops), Nsteps = nsteps, omega = 0) #6.5 hours

#################### LOOP OVER N ####################
rejection <- numeric(length(desired_indices))
for (memory in 1:n_loops) {
  #memory <- 1
  print(paste("Loop",memory, "out of ",n_loops))
  
  #Heston
  Heston <- sim.heston(settings)
  Heston_jump <- sim.addjump(Heston, burst_time = 0.5, interval_length = 0.05, c_1 = c_1, alpha = alpha)
  path <- Heston_jump
  
  #Create dy
  Y <- t(as.matrix(path$Y))
  dy <- diff(Y)
  path$Y <- t(dy)
  
  #T-estimator
  mu_hat <- est.mu.mat.2.0(data = path, hd)$mu[,desired_indices]
  sigma_hat_2 <- est.sigma.raw.mat.2.0(data = path, hd)$sig[,desired_indices]
  T_hat <- sqrt(hd/K2) * mu_hat/sqrt(sigma_hat_2)
  
  T_hat[is.na(T_hat)] <- 0 #There MIGHT be negtive sigma-hat but it is highly unlikely
  
  rejection <- rejection + apply((T_hat > threshold),2,sum)
}

rejection_perc <- rejection/Npaths*100

x_min <- -(60*10):(60*30)/60
normal_distribution <- 2.5
max_rejection <- (1-pnorm(sqrt(threshold^2-1/K2)))*100

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
