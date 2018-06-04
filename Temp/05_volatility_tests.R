#Module takes 13 seconds to run
setwd(Sys.getenv("masters-thesis"))
library(ggplot2)
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")

#seed
set.seed(100)

#Heston settings
settings <- sim.setup(mat=6.5/(24*7*52), Npath = 2, omega = 0) #6.5 hours

#Burst settings
alpha <- 0.75
beta <- 0.4
c_1 <- 0.15 #0.025% for lowest alpha. 4% for highest alpha
c_2 <- 0.016

#Get results with and without bursts. OBS: I have replaced c_1=3 with c_1 = 0.3
Heston <- sim.heston(settings)
Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = c_2, beta = beta,reverse = F, recenter = F)
Heston_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05,c_1 = c_1,alpha=alpha, reverse = F)
Heston_jump <- sim.addjump(Heston_vb,burst_time=0.5,interval_length=0.05, c_1 = c_1, alpha = alpha)


#Get a single path
path <- 2

Heston_path <- all_simulations[[1]]$Y[2,] #a = 0, b = 0, J = 0
vb_path <- all_simulations[[29]]$Y[2,] #a = 0, b = 0.4, J = 0
vbdb_path <- all_simulations[[35]]$Y[2,] #a = 0.75, b = 0.4, J = 0
jump_path <- all_simulations[[34]]$Y[2,] #a = 0.75, b = 0.4, J = 1

Heston_path <- Heston$Y[1,] #a = 0, b = 0, J = 0
vb_path <- Heston_vb$Y[1,] #a = 0, b = 0.4, J = 0
vbdb_path <- Heston_vbdb$Y[1,] #a = 0.75, b = 0.4, J = 0
jump_path <- Heston_jump$Y[1,] #a = 0.75, b = 0.4, J = 1

#Set time index for plot
time_hours <- Heston$time*(52*7*24)
time_index <- Heston$time/settings$mat

#Only plot from 0.35 to 0.65 
index <- ((time_index >=0.35) & (time_index <=0.65))
Heston_path <- Heston_path[index]
vb_path <- vb_path[index]
vbdb_path <- vbdb_path[index]
jump_path <- jump_path[index]
time_index <- time_index[index]
time_hours <- time_hours[index]


#----------------------------- NICE PLOT (NIELS' REGRESSION CODE) -----------------------------  

df_y_values <- data.frame(cbind(Heston_path,vb_path, vbdb_path, jump_path))
names <- c("Heston","+volatility burst", "+drift & volatility burst", "+jump")

tmp <- list()
for (i in 1:4) {
  tmp[[i]] <- data.frame(
    x_akse = time_hours,
    y_akse = df_y_values[,i],
    rx = names[i]
  )
}
tmp <- do.call(rbind, tmp)

ggplot(tmp, aes(x_akse, y_akse,color = rx)) +
  geom_line() +
  scale_color_manual("", values = c("black", "red", "blue", "purple")) +
  xlab("time in hours") + ylab("log-return")  + 
  #theme(legend.position = c(0.05,0),
  #      legend.justification = c(0.1,0)) +
  #      plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Log-return of asset") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))




#Alpha and beta calculations

#alpha = 0.55
c_1*(10/(60*24*7*52))^(1-0.55)/(1-0.55) #0.5%
#alpha = 0.75
c_1*(10/(60*24*7*52))^(1-0.75)/(1-0.75) #8% (extreme)

#sigma2_avg = 0.093%


