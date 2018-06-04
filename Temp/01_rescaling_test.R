setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")
source("Simulation/add_all.R")
source("Estimation/rescaling.R")

############ CHECK IF REJECTION IS IN FACT 5% ON EVERY T (WITH AND WITHOUT SCALING)

# ESTIMATION PARAMETERS
heston_params <- sim.setup()
h_mu <- 5/(52*7*24*60) #5 min as Christensen
ratio <- 15
lag <- 10
threshold <- qnorm(0.975)

# Burn a single mu in
n <- heston_params$Nsteps
mat <- heston_params$mat
dt <- mat/n
n_burn <- h_mu / dt
desired_indices <- seq(from = n_burn, to = 23400, by = 5) #Burn a mu in


#### LOOP BECAUSE OF LACK OF MEMORY, SAVE BOTH T WITHOUT AND WITH SCALING

Npaths <- 1000
n_loops <- ceiling(Npaths/200)
output_list_T <- matrix(nrow = n_loops, ncol = length(desired_indices))
output_list_T_rescaled <- output_list_T

for (memory in 1:n_loops) {
    #memory <- 1
  
    ### Keep track
    print(paste("Loop",memory, "out of ",n_loops))

    ### SIMULATE ALL PATHS
    heston_params$Npath <- ceiling(Npaths/n_loops)
    Heston <- sim.heston(heston_params)
    
    ### Calculate T's both with and without
    #dy
    Heston$Y <- t(diff(t(as.matrix(Heston$Y))))
    #T
    mu_hat <- sqrt(h_mu)*est.mu.mat.2.0(data = Heston, hd = h_mu)$mu[,desired_indices]#,t.index = t.index)$mu
    sigma_hat2 <- est.sigma.mat.3.0(data = Heston, hv = ratio*h_mu, lag = lag)$sig[,desired_indices]#,t.index = t.index,lag = lag)$sig
    T_hat <- mu_hat/sqrt(sigma_hat2)

    #Rescaling
    mu <- est.rescale.mu(mu_matrix = mu_hat, time_points = Heston$time[desired_indices], t_beginning = Heston$time[1], h_mu = h_mu)
    sigma2 <- est.rescale.sigma(sigma2_matrix = sigma_hat2, time_points = Heston$time[desired_indices], t_beginning = Heston$time[1], h_sigma = ratio*h_mu)
    T_hat_rescaled <- mu/sqrt(sigma2)
    
    #Calculate threshold (for every index)
    perc_above_T <- numeric(length = ncol(T_hat))
    perc_above_T_rescaled <- perc_above_T
    
    for (i in 1:ncol(T_hat)) {
      #i <- 1
      perc_above_T[i] <- sum((abs(T_hat[,i]) > threshold)) / nrow(T_hat)
      perc_above_T_rescaled[i] <- sum((abs(T_hat_rescaled[,i]) > threshold)) / nrow(T_hat_rescaled)
    }
  
    output_list_T[memory,] <- perc_above_T
    output_list_T_rescaled[memory,] <- perc_above_T_rescaled
}

### Take mean accross memories
T_non_scaled <- apply(output_list_T,2,mean)
T_rescaled <- apply(output_list_T_rescaled,2,mean)


### Plot
#dataframe
df_values <- c(T_non_scaled,T_rescaled,rep(0.05,length(T_non_scaled)))
df_names <- c(rep("Non-scaled", length(T_non_scaled)), rep("Re-scaled",length(T_rescaled)), rep("5%", length(T_non_scaled)))
df_time_points <- c(rep(Heston$time[desired_indices],3))
df_plot <- as.data.frame(cbind(df_values,df_names,df_time_points))
colnames(df_plot) <- c("values", "Estimator", "time_points")

#Change from factors to values
df_plot$values <- as.numeric(as.character(df_plot$values))
df_plot$time_points <- as.numeric(as.character(df_plot$time_points))*(52*7*24)

#plot
ggplot(df_plot, aes(time_points,values,color = Estimator),alpha = I(0.2)) + 
  geom_line()
  
