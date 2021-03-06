setwd(Sys.getenv("masters-thesis"))
library(latex2exp)
source("Module/table1functions.R")
source("Simulation/add_all.R")
source("Estimation/estimates_reloaded.R")


####### ESTIMATION PARAMETERS
heston_params <- sim.setup(omega = 0)
h_mu <- 5/(52*7*24*60)
threshold <- qnorm(0.975)
K2 <- 0.5

#Burn a single mu in:
n <- heston_params$Nsteps
t.index <- ceiling((n+1)/2)

# BURST SETTINGS (a = 0, 0.55, 0.65, 0.75 | b = 0.0, 0.1, 0.2, 0.3, 0.4), c_1 and c_2 from Batman
alpha <- 0.55

percentage_list <- seq(0,1, by = 0.004)
c_1_list <- -(1-alpha)*0.01*percentage_list/(10/(60*24*7*52))^(1-alpha) #-c_1 creates positive jump (bad coded)

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
Npaths <- 1000 #Takes approx. a second per path (because it has to estimate T for 35 processes w. 3 different bandwidths)
n_loops <- ceiling(Npaths/50) #After 50 it just scales linearly if not slower
output_list <- list()
set.seed(100)
N <- ceiling(Npaths/n_loops)
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
      print(process)
      
      #Calculate dy
      all_simulations[[process]]$Y <- t(diff(t(as.matrix(all_simulations[[process]]$Y))))
      #h <- all_simulations[[process]]$Y
      #Calculate T
      mu_hat <- est.mu.mat.2.0(data = all_simulations[[process]], hd = h_mu)$mu[,t.index]
      sigma2_hat <- est.sigma.raw.mat.2.0(data = all_simulations[[process]], hv = h_mu)$sig[,t.index]
      T_hat <- sqrt(h_mu/K2) * mu_hat/ sqrt(sigma2_hat)

        
      #Save percentage rejection in list
      temp[memory,process] <- mean(T_hat > threshold)
      }
}
print(Sys.time()-p0)

### Take mean accross memories (assuming same number of paths in every memory loop)
results_mean <- apply(na.omit(temp),2,mean)*100 #na.omit because i stop it before it is finished

#plot
normal_distribution <- 2.5
max_rejection <- (1-pnorm(sqrt(q^2-c^2)))*100

#plot1
x <- rep(percentage_list,3)

y_list <- c(results_mean,
            rep(normal_distribution,length(c_1_list)),
            rep(max_rejection,length(c_1_list)))

names <- c(rep('rejection',length(results_mean)),
           rep('2.5 %',length(results_mean)),
           rep('theoretical max', length(results_mean)))


plot_df <- data.frame(do.call('cbind',list(x,y_list,names)))


ggplot(plot_df,aes(x,y_list, color = names)) +
  geom_line() +
  xlab("Jump size in %") +
  ylab(TeX('$ P(T > q_{0.975})$')) +
  ggtitle(TeX('$ P(T > q_{0.975})$ for various jump sizes')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  theme(legend.title = element_blank())

#save(plot_df, file="Figures2/Saved_data_for_plots/04_extra_diff_jumps.Rda")



#plot1
perc_list <- percentage_list[1:round((length(percentage_list)+1)/2,0)] 

x <- rep(perc_list,3)

y_list <- c(results_mean[1:length(perc_list)],
            rep(normal_distribution,length(perc_list)),
            rep(max_rejection,length(perc_list)))

names <- c(rep('rejection',length(perc_list)),
           rep('2.5 %',length(perc_list)),
           rep('theoretical max', length(perc_list)))


plot_df <- data.frame(do.call('cbind',list(x,y_list,names)))


ggplot(plot_df,aes(x,y_list, color = names)) +
  geom_line() +
  xlab("Jump size in %") +
  ylab(TeX('$ P(T > q_{0.975})$')) +
  ggtitle(TeX('$ P(T > q_{0.975})$ for various jump sizes')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  theme(legend.title = element_blank())

