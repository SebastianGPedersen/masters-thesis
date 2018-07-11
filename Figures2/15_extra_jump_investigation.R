library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")
source("Simulation/add_all.R")


####### ESTIMATION PARAMETERS
heston_params <- sim.setup()
h_mu <- 5/(52*7*24*60)
ratio <- 12
lag <- 10

#Burn a single minute in:
n <- heston_params$Nsteps
mat <- heston_params$mat
dt <- mat/n
one_min <- 1/(52*7*24*60)
n_burn <- one_min / dt
t.index <- seq(from = n_burn, to = 23400, by = 5) #Burn a volatility bandwidth (note 10 in Christensen)

# BURST SETTINGS (a = 0, 0.55, 0.65, 0.75 | b = 0.0, 0.1, 0.2, 0.3, 0.4), c_1 and c_2 from Batman
alpha <- 0.55

percentage_list <- seq(0,1, by = 0.01)
c_1_list <- (1-alpha)*0.01*percentage_list/(10/(60*24*7*52))^(1-alpha)

#Get burst settings as a list (uses sim.burstsetting for standardization)
burstsettings <- list(length(c_1_list))
for (c_1 in 1:length(c_1_list)) {
  #beta_index <- 1
  burstsettings[[c_1]] <- sim.burstsetting(jump = T,
                                           alpha = alpha,
                                           c_1 = c_1_list[c_1],
                                           beta = 0,
                                           c_2 = 0)
}

#Evt. reverse
#### LOOP BECAUSE OF LACK OF MEMORY
Npaths <- 500 #Takes approx. a second per path (because it has to estimate T for 35 processes w. 3 different bandwidths)
n_loops <- ceiling(Npaths/50) #After 50 it just scales linearly if not slower
output_list <- list()
set.seed(100)
N <- ceiling(Npaths/n_loops)
Tstar <- numeric(N);  rho<-numeric(N);
temp <- matrix(nrow = n_loops,ncol = length(burstsettings))

p0 <- Sys.time()
for (memory in 1:n_loops) {
    #memory <- 1
  
    ### Keep track
    print(paste("Loop",memory, "out of ",n_loops))
    print(paste("Expected time left:", round(Npaths-(memory-1)*Npaths/n_loops,0),"seconds")) #One second per path
    
    ### SIMULATE ALL PATHS
    heston_params$Npath <- N
    Heston <- sim.heston(heston_params)
    
    all_simulations <- sim.add_all(Heston = Heston, burst_args = burstsettings)
    
    for (process in 1:length(all_simulations)) {
      #process <- 1
      
      #Calculate dy
      all_simulations[[process]]$Y <- t(diff(t(as.matrix(all_simulations[[process]]$Y))))
      
      #Calculate T
      mu_hat <- est.mu.mat.2.0(data = all_simulations[[process]], hd = h_mu, bandwidth_rescale = T)$mu[,t.index]
      sigma2_hat <- est.sigma.mat.3.0(data = all_simulations[[process]], hv = ratio*h_mu, bandwidth_rescale = T)$sig[,t.index]
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
      temp[memory,process] <- mean(Tstar>=z)
      }
}
print(Sys.time()-p0)

### Take mean accross memories (assuming same number of paths in every memory loop)
results_mean <- apply(temp,2,mean)*100

#plot
plot_df <- data.frame(percentage = percentage_list,
                      rejection = results_mean,
                      color = rep("Jumps",length(percentage_list)))

five_perc_df <- data.frame(percentage = percentage_list,
                           rejection = rep(5,length(percentage_list)),
                           color = rep("5%",length(percentage_list)))

plot_df <- rbind(plot_df,five_perc_df)

ggplot(plot_df, aes(percentage,rejection,color = color)) +
  geom_line() +
  xlab("Jump size in %") +
  ylab(TeX('$ P(T^*> q_{0.95})$')) +
  ggtitle(TeX('$ P(T^*> q_{0.95})$ for various jump sizes')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  theme(legend.title = element_blank())

