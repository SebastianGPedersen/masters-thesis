#Module takes 13 seconds to run
setwd(Sys.getenv("masters-thesis"))
library(ggplot2)
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")

#seed
set.seed(100)

#Heston settings
settings <- sim.setup(mat=6.5/(24*7*52), Npath = ceiling(Npaths/n_loops), Nsteps = nsteps, omega = 0) #6.5 hours
dt <- settings$mat /nsteps

#Burst settings
q <- threshold
c <- sqrt(2)
sigma <- sqrt(settings$theta)
k_max <- c*sigma/sqrt(q^2-c^2)
J <- -k_max*sqrt(h_mu) #0.06%
alpha <- 0.8
c_1 <- (1-alpha)*J/(10/(60*24*7*52))^(1-alpha)


#Get results with and without bursts. OBS: I have replaced c_1=3 with c_1 = 0.3
BS <- sim.BlackScholes(mean = 0, sd = sigma, omega = 0, Nsteps = nsteps, Npath = 2)
BS_jump <- sim.addjump(BS, burst_time = 0.5, interval_length = 0.05, c_1 = c_1, alpha = alpha)

#Get a single path
path <- 1

Heston_path <- sim.path(path,BS)$Y
jump_path <- sim.path(path, BS_jump)$Y

#Set time index for plot
time_hours <- BS$time *(52*7*24)
time_index <- BS$time/settings$mat

#Only plot from 0.35 to 0.65 
index <- ((time_index >=0.35) & (time_index <=0.65))
Heston_path <- Heston_path[index]
jump_path <- jump_path[index]
time_index <- time_index[index]
time_hours <- time_hours[index]


#----------------------------- NICE PLOT (NIELS' REGRESSION CODE) -----------------------------  

df_y_values <- data.frame(cbind(Heston_path, jump_path))
names <- c("Black Scholes", "+jump")

tmp <- list()
for (i in 1:2) {
  tmp[[i]] <- data.frame(
    x_akse = time_hours,
    y_akse = df_y_values[,i],
    rx = names[i]
  )
}
tmp <- do.call(rbind, tmp)

ggplot(tmp, aes(x_akse, y_akse,color = rx)) +
  geom_line() +
  scale_color_manual("", values = c("black", "red")) +
  xlab("time in hours") + ylab("log-return")  + 
  #theme(legend.position = c(1,0),
  #      legend.justification = c(1,0),
  #      plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Log-return of asset") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))



dx <- jump_path[2:length(jump_path)]-jump_path[1:(length(jump_path)-1)]
sd(dx)