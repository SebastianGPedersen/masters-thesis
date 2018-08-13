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
omega2 <- 2.64*10^(-10) #What Mathias wrote
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400 /7
mat <- 6.5/(24*7*52) /7 #I cheat a little because it takes too long otherwise. Since it is only evaluated in a single point, it doesn't matter as long as there are enough bandwidths of obs.
dt <- mat/n
Npaths <- 1000
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
#noise_ratio <- omega/(sigma*sqrt(dt))

#Because of lack of memory, it is done in loops
n_loops <- ceiling(Npaths/100)

#List to final values
hd_list <- seq(1,20,by = 1)/ (2*60*24*7*52)
var_bias_temp <- matrix(nrow = n_loops,ncol = length(hd_list))

for (memory in 1:n_loops) {
  #memory <- 1
  print(memory)
  temp_paths <- Npaths / n_loops
  BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = omega, Nsteps = n, Npath = temp_paths)

  #Create dy
  BS$Y <- t(diff(t(as.matrix(BS$Y))))
  
  #################### LOOP OVER H_D ####################
  
  #H_d from 1min to 15min
  for (hd in 1:length(hd_list)) {
    
    #hd <- 2
    desired_index <- n-1 #Takes last index, so K uses as many points as possible
    
    mu_hat <- sqrt(hd_list[hd])/sqrt(K2*sigma^2)*est.mu.mat(data = BS, hd = hd_list[hd], t.index = desired_index)$mu
    var_bias_temp[memory,hd] <- mean(mu_hat^2) #We know it has mean zero, so E(mu^2) is more precise than Var(mu)
  }
}

#The total variance is just the mean of the variances (because the loops have same size)
var_bias <- numeric(length = length(hd_list))
for (hd in 1:length(hd_list)) {
  var_bias[hd] <- mean(var_bias_temp[,hd])
}

#Compute the the bias from the last noise term
noise <- 1+omega^2 / (hd_list*K2*sigma2)


#################### PLOT ####################
hd_minutes <- hd_list*(60*24*7*52)
data <- data.frame(hd = hd_minutes, target = (1:length(hd_minutes)*0)+1, var_bias = var_bias, var_bias_noise = noise)

ggplot() +
  geom_line(data=data, aes(x=hd, y=var_bias, col = "Var_mu"), size = 1) +
  geom_line(data=data, aes(x=hd, y=var_bias_noise, col = "var-noise"), size = 1) +
  geom_line(data=data, aes(x=hd, y=target, col = "1"), size = 1) +
  xlab(TeX("Bandwidth in minutes")) + ylab("Value") +
  ggtitle('Bias of drift estimator') +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  scale_color_discrete(name = "Expression",
                       labels = unname(TeX(
    c("1",'Var $\\left( \\frac{\\sqrt{h_n} \\mu_t^n}{\\sqrt{K_2 \\sigma^2}} \\right)$',
      "$1+\\frac{\\omega^2}{h_n K_2 \\sigma^2$"))))

