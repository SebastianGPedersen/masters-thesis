setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")
source("Simulation/add_all.R")


####### ESTIMATION PARAMETERS
heston_params <- sim.setup()
h_mu <- 5/(52*7*24*60)
ratio <- 12
lag <- 10

ten_minutes <- 10 / (6.5*60) #Starts immeiately
two_point_five_min <- 2.5 / (6.5*60) #Starts after 3hours 5min

burst_times <- seq(ten_minutes,0.5,by = two_point_five_min)
begin_times <- seq(0,0.5-ten_minutes,by = two_point_five_min)


#Burn a single mu in:
n <- heston_params$Nsteps
mat <- heston_params$mat
dt <- mat/n
n_burn <- h_mu / dt
t.index <- seq(from = n_burn, to = 23400, by = 5) #Burn a volatility bandwidth (note 10 in Christensen)

# BURST SETTINGS (a = 0, 0.55, 0.65, 0.75 | b = 0.0, 0.1, 0.2, 0.3, 0.4), c_1 and c_2 from Batman
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
    c_2 <- sqrt((1-2*beta)*(0.00093*0.25*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.2) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.5*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.3) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.75*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.4) {
    c_2 <- sqrt((1-2*beta)*(0.00093*1*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  return(c_2)
}


#Create (vol, drift, jump)*burst_times

burstsettings <- lapply(burst_times, function(x) sim.burstsetting(burst_time = x,jump = F,alpha = 0, c_1 = 0, beta = 0.3, c_2 = c_2_func(0.3)))
burstsettings <- append(burstsettings, 
                        lapply(burst_times, function(x) sim.burstsetting(burst_time = x,jump = F,alpha = 0.55, c_1 = c_1_func(0.55), beta = 0, c_2 = 0)))
burstsettings <- append(burstsettings, 
                        lapply(burst_times, function(x) sim.burstsetting(burst_time = x,jump = T,alpha = 0.55, c_1 = c_1_func(0.55), beta = 0, c_2 = 0)))

#225 bursts processes

#Evt. reverse
#### LOOP BECAUSE OF LACK OF MEMORY
Npaths <- 200 #Takes approx. a second per path (because it has to estimate T for 35 processes w. 3 different bandwidths)
n_loops <- ceiling(Npaths/50) #After 50 it just scales linearly if not slower
output_list <- list()
set.seed(100)
N <- ceiling(Npaths/n_loops)
Tstar <- numeric(N);  rho<-numeric(N);

p0 <- Sys.time()
results <- matrix(0,nrow = 3,ncol = length(burst_times))

for (memory in 1:n_loops) {
    #memory <- 1
  
    ### Keep track
    print(paste("Loop",memory, "out of ",n_loops))
    print(paste("Expected time left:", 3*round(Npaths-(memory-1)*Npaths/n_loops,0),"seconds")) #One second per path
    
    ### SIMULATE ALL PATHS
    heston_params$Npath <- N
    Heston <- sim.heston(heston_params)
    
    all_simulations <- sim.add_all2(Heston = Heston, burst_args = burstsettings)
    
    for (process in 1:3) {
      for (burst_time in 1:length(burst_times)) {
      #process <- 1
      #burst_time <- 1
      print(paste(memory,(process-1)*length(burst_times)+burst_time))
      path <- all_simulations[[(process-1)*length(burst_times)+burst_time]]
      #Calculate dy
      path$Y <- t(diff(t(as.matrix(path$Y))))
      
      #Calculate T
      mu_hat <- est.mu.mat.2.0(data = path, hd = h_mu, bandwidth_rescale = T)$mu[,t.index]
      sigma2_hat <- est.sigma.mat.3.0(data = path, hv = ratio*h_mu, bandwidth_rescale = T)$sig[,t.index]
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
      results[process,burst_time] <- results[process,burst_time] + mean(Tstar>=z)
      }
    }
}
print(Sys.time()-p0)

### Take mean accross memories (assuming same number of paths in every memory loop)
results_mean <- results / n_loops


### Restructure
names <- c("vol burst", "drift burst", "jump")
times_in_hours <- begin_times*6.5

#Create a single data_frame
plot_data_frame <- data.frame(x_axis = times_in_hours,
                              rejection = rep(0.05,length(burst_times)),
                              process = "5%")

for (i in 1:nrow(results)){
  new_data_frame <- data.frame(x_axis = times_in_hours,
                               rejection = results_mean[i,],
                               process = names[i])
  plot_data_frame <- rbind(plot_data_frame,new_data_frame)
}


##### PLOT #####
ggplot(plot_data_frame, aes(x_axis, rejection, color = process)) +
  geom_line() +
  xlab("Beginning of burst in hours") + ylab('Rejection percentage') +
  ggtitle("Rejection of T-estimator at different burst times") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))



