library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/add_all.R")
source("Estimation/estimates_reloaded.R")
source("Estimation/estimates_revolution.R")
source("Estimation/rho.R")
source("estimation/teststat.R")


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
one_min <- 5 / (52*7*24*60)
n_burn <- one_min / dt
t.index <- seq(from = n_burn, to = n, by = 5) #Burn a volatility bandwidth (note 10 in Christensen)

#Burst settings1
alpha_list <- c(0.55,0.65,0.75)

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


burstsettings <- list(length(alpha_list))
for (alpha in 1:length(alpha_list)) {
  burstsettings[[alpha]] <- sim.burstsetting(jump = T,
                                             alpha = alpha_list[alpha],
                                             c_1 = c_1_func(alpha_list[alpha]),
                                             beta = 0,
                                             c_2 = 0)
}


#The ratio parameters
ratio_list <- seq(1,30,by = 1) #From 1 to 20
h_mu_list <- c(2,5,10) / (52*7*24*60)


#Because of lack of memory, it is done in loops
Npaths <- 300 #15min
n_loops <- ceiling(Npaths/50) #After 50 it just scales linearly if not slower
output_list <- list()
N <- ceiling(Npaths/n_loops)

#Initialize lists used later
Tstar <- numeric(N);  rho<-numeric(N);
temp <- matrix(nrow = length(burstsettings), ncol = length(ratio_list))
my_list <- rep(list(temp),n_loops)
output_list <- rep(list(my_list),length(h_mu_list))

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
  
  #For all bandwidths, save the rejection avg. for every jump size
  for (h_mu in 1:length(h_mu_list)){
    #h_mu <- 1
    for (process in 1:length(all_simulations)) {
      #process <- 1
      mu_hat <- est.mu.mat.2.0(data = all_simulations[[process]], hd = h_mu_list[h_mu], bandwidth_rescale = T)$mu[,t.index]
      
      for (ratio in 1:length(ratio_list)) {
        #ratio <- 5
          
        #Estimate T
        sigma2_hat <- est.sigma.mat.3.0(data = all_simulations[[process]], hv = ratio_list[ratio]*h_mu_list[h_mu], bandwidth_rescale = T)$sig[,t.index]
        T_hat <- sqrt(h_mu_list[h_mu]) * mu_hat/ sqrt(sigma2_hat)
        
        #Calculate T_max and thresholds
        for(subpath in 1:N){
          Tstar[subpath] <- max(abs(T_hat[subpath,]))
          # fit rho
          rho[subpath] <- est.rho.MLE(T_hat[subpath,])
        }
  
        m <- rep(ncol(T_hat),N)
        #Save percentage rejection in list
        temp[process,ratio] <- mean(Tstar>=z)
      }
    }
    output_list[[h_mu]][[memory]] <- temp
  }
}
print(Sys.time()-p0)

###Take mean over memory loops
temp_matrix <- matrix(0,nrow = length(all_simulations), ncol = length(ratio_list))
output_mean <- rep(list(temp_matrix),length(h_mu_list))

for (h_mu in 1:length(h_mu_list)) {
  for (memory in 1:n_loops) {
    output_mean[[h_mu]] <- output_mean[[h_mu]] + output_list[[h_mu]][[memory]]
  }
  output_mean[[h_mu]] <- output_mean[[h_mu]] / n_loops*100
}

### Find upper and lower for every h_mu
for (h_mu in 1:length(h_mu_list)) {
  output_mean[[h_mu]] <- apply(output_mean[[h_mu]],2,sort)
}



### Restructure
h_mins <- round(h_mu_list*(52*7*24*60),0)
  
#Create a single data_frame
plot_data_frame <- data.frame(ratios = ratio_list ,
                              lower= rep(5,length(ratio_list)),
                              mean= rep(5,length(ratio_list)),
                              upper= rep(5,length(ratio_list)),
                              mu_bandwidth = "5%")

for (i in 1:length(output_mean)){
  new_data_frame <- data.frame(ratios = ratio_list,
                               lower = output_mean[[i]][1,],
                               mean = output_mean[[i]][2,],
                               upper= output_mean[[i]][3,],
                               mu_bandwidth = paste("h_mu = ",h_mins[i]))
  plot_data_frame <- rbind(plot_data_frame,new_data_frame)
}


##### PLOT #####
ggplot(plot_data_frame, aes(ratios, mean, color = mu_bandwidth)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = mu_bandwidth), alpha = 0.3) +
  xlab("Bandwidth ratio") + ylab('Rejection percentage') +
  ggtitle("Rejection of T-estimator with jumps") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
  
save(plot_data_frame, file="Figures2/Saved_data_for_plots/13_bandwidth_ratio.Rda")

