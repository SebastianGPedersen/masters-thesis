library(ggplot2)
library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/BlackScholes.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")
source("Estimation/estimates_reloaded.R")

set.seed(100)

#################### PARAMETERS THAT DON'T CHANGE ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400 /7
mat <- 6.5/(24*7*52) /7
dt <- mat/n
Npaths <- 1000
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sqrt(dt)*sigma)) #close to 1/2

#Because of lack of memory, it is done in loops
n_loops <- ceiling(Npaths/100)

#Initialize lists to contain final values
hd_list <- seq(1,20,by = 1)/ (2*60*24*7*52) #From 30sec to 10min?
var_mu_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))
sig_temp1 <- matrix(nrow = n_loops,ncol = length(hd_list))
var_T_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))

for (memory in 1:n_loops) {
  #memory <- 1
  temp_paths <- Npaths / n_loops
  BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = omega, Nsteps = n, Npath = temp_paths)
  
  #### Calculate dY ####
  BS$Y <- t(diff(t(as.matrix(BS$Y))))

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
    lag <- 10
    
    ### Calculate sigma, mu and t ###
    desired_index <- n-1 #Takes last index, so K uses as many points as possible
    
    #With matrices
    sig_hat1 <- 1/(K2*sigma^2) * est.sigma.mat(BS, hv = hd_list[hd], t.index = desired_index, lag = lag)$sig
    
    mu_hat <- sqrt(hd_list[hd])/sqrt(K2*sigma^2)*est.mu.mat(data = BS, hd = hd_list[hd], t.index = desired_index)$mu
    T_hat <- mu_hat/sqrt(sig_hat1)
    
    ### Save results in matrices
    var_mu_bias_temp[memory,hd] <- mean(mu_hat^2) #We know true mean is zero so uses squared mean instead of var
    sig_temp1[memory,hd] <- mean(sig_hat1)

    var_T_bias_temp[memory,hd] <- var(na.omit(T_hat)) #Mean is not zero, so uses 'var'.
  }
}

#################### TAKE MEAN OVER LOOPS  ####################
var_mu_bias <- numeric(length = length(hd_list))
sig_bias1 <- var_mu_bias
var_T_bias <- var_mu_bias

for (hd in 1:length(hd_list)) {
  var_mu_bias[hd] <- mean(var_mu_bias_temp[,hd])
  sig_bias1[hd] <- mean(sig_temp1[,hd])
  
  var_T_bias[hd] <- mean(var_T_bias_temp[,hd])
}

noise <- 1/hd_list * omega^2 / (K2*sigma^2) + 1

#################### PLOT ####################
hd_minutes <- hd_list*(60*24*7*52)
plot_data_frame <- data.frame(hd = hd_minutes, 
                              target = (1:length(hd_minutes)*0)+1, 
                              mu_bias = var_mu_bias, 
                              sigma_bias1 = sig_bias1, 
                              T_bias = var_T_bias, 
                              last_noise_bias = noise)
#load("Figures2/Saved_data_for_plots/09_bias_of_all1.Rda")

ggplot() +
  geom_line(data=plot_data_frame, aes(x=hd, y=mu_bias, col = "mu_bias"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=last_noise_bias, col = "Last_noise_bias"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=target, col = "1"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=sigma_bias1, col = "sigma_bias"), size = 1) +
  geom_line(data=plot_data_frame, aes(x=hd, y=T_bias, col = "T_bias"), size = 1) +
  xlab(TeX("Bandwidth in minutes")) + ylab("Value") +
  ggtitle(TeX('Bias of volatility and $T$-estimator')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  scale_color_discrete(name = "Expression",
                       labels = unname(TeX(
                         c("1","$1+\\frac{\\omega^2}{h_n K_2 \\sigma^2$","$T_t^n$",
                           'Var $\\left( \\frac{\\sqrt{h_n} \\mu_t^n}{\\sqrt{K_2 \\sigma^2}} \\right)$',
                           "$E \\left( \\frac{\\hat{\\Sigma}_t^n}{K_2 \\sigma_{t-}^2} \\right)$"))))
#Save dataframe for later
#save(plot_data_frame, file="Figures2/Saved_data_for_plots/09_bias_of_all1.Rda")
