library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/BlackScholes.R")
source("Estimation/estimates_reloaded.R")
source("Estimation/estimates_revolution.R")
source("Estimation/laglength.R")

p0 <- Sys.time()
#################### PARAMETERS THAT DON'T CHANGE ####################
omega2 <- 2.64*10^(-10)*25
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400  /7
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 1000
sigma2 <- 0.0457/25
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sqrt(dt)*sigma)) #close to 1/2

#Because of lack of memory, it is done in loops
n_loops <- ceiling(Npaths/200) 
N <- ceiling(Npaths/n_loops)

#Initialize lists to contain final values
hd_list <- seq(1,20,by = 1)/ (2*60*24*7*52)
var_mu_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))
sig_temp1 <- matrix(nrow = n_loops,ncol = length(hd_list))
sig_temp2 <- matrix(nrow = n_loops,ncol = length(hd_list))
sig_temp3 <- matrix(nrow = n_loops,ncol = length(hd_list))
var_T_bias_temp1 <- matrix(nrow = n_loops,ncol = length(hd_list))
var_T_bias_temp2 <- matrix(nrow = n_loops,ncol = length(hd_list))
var_T_bias_temp3 <- matrix(nrow = n_loops,ncol = length(hd_list))


for (memory in 1:n_loops) {
  set.seed(1000*memory)
  
  #memory <- 1
  BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = omega, Nsteps = n, Npath = N)
  
  #### Calculate dY ####
  dY <- t(diff(t(BS$Y))) #Has to be transposed for 'diff' to work as desired
  data <- list(Y = dY, time = BS$time)

  #################### LOOP OVER BANDWIDTH ####################
  for (hd in 1:length(hd_list)) {
    print(paste0("memory = ",memory, ", hd = ",hd,sep = ""))
    #hd <- 1

    ### Lag selection ###
    #With Cristensen et. al. code
    #lag <- laglength(dY[1,],hd_list[hd]/dt) #The second is from the fact that: sum(K) -> h_n/dt * int(0,1,K(x))
    #print(paste("lag = ",lag,sep = ""))
    
    #With eq. (3) in Barndorff-Nielsen (2009)
    #c_star <- 3.5134
    #xi <- sqrt(omega^2/sqrt(mat^2*sigma^4)) #If mat increases.
    #lag <- c_star*xi^(4/5)*(n)^(3/5)
    
    #With my heurestic guess from my proof
    #lag <- c_star*xi^(4/5)*(hd_list[hd]/dt)^(1.1))  #It actually seems to work quide well

    #With constant lag
    lags <- c(5,10,100)
    
    ### Calculate sigma, mu and t ###
    desired_index <- n-1 #Takes last index, so K uses as many points as possible
    
    #With matrices
    sig_hat1 <- 1/(K2*sigma^2) * est.sigma.mat.2.0(data, hv = hd_list[hd], lag = lags[1])$sig[,desired_index]
    sig_hat2 <- 1/(K2*sigma^2) * est.sigma.mat.3.0(data, hv = hd_list[hd], lag = lags[2])$sig[,desired_index]
    sig_hat3 <- 1/(K2*sigma^2) * est.sigma.mat.3.0(data, hv = hd_list[hd], lag = lags[3])$sig[,desired_index]
    
    mu_hat <- sqrt(hd_list[hd])/sqrt(K2*sigma^2)*est.mu.mat.2.0(data = data, hd = hd_list[hd])$mu[,desired_index]
    
    T_hat1 <- mu_hat/sqrt(sig_hat1)
    T_hat2 <- mu_hat/sqrt(sig_hat2)
    T_hat3 <- mu_hat/sqrt(sig_hat3)
    
    ### Save results in matrices
    var_mu_bias_temp[memory,hd] <- mean(mu_hat^2) #We know true mean is zero so uses squared mean instead of var
    sig_temp1[memory,hd] <- mean(sig_hat1)
    sig_temp2[memory,hd] <- mean(sig_hat2)
    sig_temp3[memory,hd] <- mean(sig_hat3)
    
    var_T_bias_temp1[memory,hd] <- var(na.omit(T_hat1)) #Mean is not zero, so uses 'var'.
    var_T_bias_temp2[memory,hd] <- var(na.omit(T_hat2)) #Mean is not zero, so uses 'var'.
    var_T_bias_temp3[memory,hd] <- var(na.omit(T_hat3)) #Mean is not zero, so uses 'var'.
    
  }
}

#################### TAKE MEAN OVER LOOPS  ####################
var_mu_bias <- numeric(length = length(hd_list))
sig_bias1 <- var_mu_bias
sig_bias2 <- var_mu_bias
sig_bias3 <- var_mu_bias
var_T_bias1 <- var_mu_bias
var_T_bias2 <- var_mu_bias
var_T_bias3 <- var_mu_bias

for (hd in 1:length(hd_list)) {
  var_mu_bias[hd] <- mean(var_mu_bias_temp[,hd])
  sig_bias1[hd] <- mean(sig_temp1[,hd])
  sig_bias2[hd] <- mean(sig_temp2[,hd])
  sig_bias3[hd] <- mean(sig_temp3[,hd])
  
  var_T_bias1[hd] <- mean(var_T_bias_temp1[,hd])
  var_T_bias2[hd] <- mean(var_T_bias_temp2[,hd])
  var_T_bias3[hd] <- mean(var_T_bias_temp3[,hd])
  
}


#################### PLOT ####################
noise <- 1/hd_list * omega^2 / (K2*sigma^2) + 1

hd_minutes <- hd_list*(60*24*7*52)
plot_data_frame <- data.frame(hd = hd_minutes, 
                              target = (1:length(hd_minutes)*0)+1, 
                              mu_bias = var_mu_bias, 
                              sigma_bias1 = sig_bias1, 
                              sigma_bias2 = sig_bias2, 
                              sigma_bias3 = sig_bias3, 
                              T_bias1 = var_T_bias1,
                              T_bias2 = var_T_bias2,
                              T_bias3 = var_T_bias3,
                              last_noise_bias = noise)

#load("Figures2/Saved_data_for_plots/09_bias_of_all1.Rda")
ggplot() +
  geom_line(data=plot_data_frame, aes(x=hd, y=mu_bias, col = "mu_bias"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=target, col = "1"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=sigma_bias1, col = "sigma_bias, lag =  5"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=sigma_bias2, col = "sigma_bias, lag = 10"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=sigma_bias3, col = "sigma_bias, lag = 100"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=T_bias1, col = "T_bias, lag =  5"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=T_bias2, col = "T_bias, lag = 10"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=T_bias3, col = "T_bias, lag = 100"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=last_noise_bias, col = "Last_noise_bias"), size = 1) +
  xlab("Bandwidth, hd") + ylab("Normalized Variance")

#Save dataframe for later
#save(plot_data_frame, file="Figures2/Saved_data_for_plots/09_bias_of_all1.Rda")

print(Sys.time()-p0)
