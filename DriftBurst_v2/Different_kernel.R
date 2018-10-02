# Function from 0 to 1 to insert in pre-avg.
parzen_kernel_func <- function(x) {
  new_x <- abs((x-1/2)*2)
  
  if (x < 1/2) {parz <- 1-6*new_x^2+6*new_x^3}
  else {parz <- 2*(1-new_x)^3}
  
  return(parz)

}


#Kommentar fra Sebastian:
#   Grunden til de får at drift er normalfordelt er fordi de sætter k_n = 3,
#   så gælder det netop at leddet psi_1 * k_n = 1 og vores skalering er irrelevant.
#   Jeg viser dette nedenstående ved at lave QQ-plot med forskellige k_n værdier


set.seed(100)

library(ggplot2)
library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/BlackScholes.R")
source("Estimation/estimates.R")
source("Estimation/pre-average.R")

#### Parameters ####
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400 
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 200
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
h_mu <- 5 / (60*24*7*52)
lag <- 8

### k_n's
k_n_list <- c(2,3,4,5,10)


#Because of lack of memory, it is done in loops
n_loops <- ceiling(Npaths/100)
desired_index <- n-1


#List to final values
T_estimator <- matrix(nrow = Npaths,ncol = length(k_n_list))


for (memory in 1:n_loops) {
  #memory <- 1
  temp_paths <- Npaths / n_loops
  
  ### Simulate BS
  BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = omega, Nsteps = n, Npath = temp_paths)
  Y <- t(as.matrix(BS$Y))
  dy <- diff(Y)  
  
  for (k_n in 1:length(k_n_list)) {
    #k_n <- 1
    print(paste0("memory = ",memory, ", k_n = ",k_n,sep = ""))
    
    ###Pre-averaging doesn't take vector
    pre_y <- matrix(NA, nrow = temp_paths, ncol = n-1) #Does NOT include 0 because BS$Y does not
    
    
    for (i in 1:temp_paths){
      pre_y[i,] <- c(rep(0,k_n_list[k_n]-2),est.NewPreAverage(dy[,i],k_n_list[k_n], kernel_func = parzen_kernel_func))
    }
    
    BS$Y <- pre_y
    
    ## Calculate mu
    mu_hat <- est.mu.mat(BS, hd = h_mu, t.index = desired_index)$mu
    sig_hat <- est.sigma.mat(BS, hv = 5*h_mu, t.index = desired_index, lag = lag)$sig
    
    #T_estimator[((memory-1)*temp_paths+1):(memory*temp_paths),k_n] <- sqrt(h_mu) * mu_hat / sqrt(K2*sig_hat)
    T_estimator[((memory-1)*temp_paths+1):(memory*temp_paths),k_n] <- sqrt(h_mu) * mu_hat / sqrt(sig_hat)
    
  }
  
}


#Re-shape to data-frame
plot_data_frame <- data.frame(do.call("rbind",
                                      list(cbind(T_estimator[,1],rep("2",dim(T_estimator)[1])), 
                                           cbind(T_estimator[,2],rep("3",dim(T_estimator)[1])),
                                           cbind(T_estimator[,3],rep("4",dim(T_estimator)[1])),
                                           cbind(T_estimator[,4],rep("5",dim(T_estimator)[1])),
                                           cbind(T_estimator[,5],rep("10",dim(T_estimator)[1])))))

colnames(plot_data_frame) <- c("T_estimator", "kn")

#Numeric and factor
plot_data_frame$T_estimator <- as.numeric(as.character(plot_data_frame$T_estimator))
plot_data_frame$kn <- as.factor(plot_data_frame$kn)


ggplot(plot_data_frame) + 
  stat_qq(aes(sample = T_estimator, colour = kn)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Quantile of standard normal") + ylab("Quantile of T") +
  ggtitle(TeX("T-estimator from Drift Burst v2 (without $K_2$)")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) +
  scale_color_discrete(name = unname(TeX("$k_n")))

save(plot_data_frame, file = "DriftBurst_v2/Saved_data_for_plots/Different_kernel.Rda")

