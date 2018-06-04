setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")
source("Simulation/add_all.R")
source("Estimation/gumbel.R")
source("estimation/rho.R")
source("estimation/teststat.R")


############ CHECK IF REJECTION IS IN FACT 5% ON EVERY T (WITH AND WITHOUT SCALING)

# ESTIMATION PARAMETERS
heston_params <- sim.setup()
h_mu <- 5/(52*7*24*60) # 5 min as Christensen
ratio <- 15 #As Christensen
lag <- 10

# Burn a single mu in
n <- heston_params$Nsteps
mat <- heston_params$mat
dt <- mat/n
n_burn <- h_mu*ratio / dt
desired_indices <- seq(from = n_burn, to = 23400, by = 60) #Burn a sigma in
m <- length(desired_indices)

#Threshold
#threshold <- qnorm(0.975)
threshold <- q95(m)

#### LOOP BECAUSE OF LACK OF MEMORY, SAVE BOTH T WITHOUT AND WITH SCALING

Npaths <- 1000
n_loops <- ceiling(Npaths/200)
output_list_T <- vector(length = n_loops)
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
    mu <- est.rescale.mu(mu_vector = mu_hat, time_points = Heston$time[desired_indices], t_beginning = Heston$time[1], bandwidth = h_mu)
    sigma2 <- est.rescale.sigma(sigma2_vector = sigma_hat2, time_points = Heston$time[desired_indices], t_beginning = Heston$time[1], bandwidth = ratio*h_mu)
    T_rescaled <- mu/sqrt(sigma2)
    
    #thresholds w. rho
    rho_vector <- apply(T_hat,1,est.rho.MLE) #Get all rho's. Seems fine
    m_vector <- rep(m,length(rho_vector)) #Repeat m's
    thresholds <- est.z_quantile(m_vector, rho_vector, 0.95)$qZm #Thresholds.
    
    #Calculate max
    T_star <- apply(abs(T_hat),1,max)
    T_star_rescaled <- apply(abs(T_rescaled),1,max)
    
    output_list_T[memory] <- sum((T_star > thresholds)) / length(T_star)
    output_list_T_rescaled[memory] <- sum((T_star_rescaled > thresholds)) / length(T_star_rescaled)
}

### Take mean accross memories
T_non_scaled <- mean(output_list_T)
T_rescaled <- mean(output_list_T_rescaled)

print(T_non_scaled)
print(T_rescaled)
