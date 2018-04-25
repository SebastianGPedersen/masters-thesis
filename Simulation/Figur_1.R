#Module takes 13 seconds to run
setwd(Sys.getenv("masters-thesis"))
library(ggplot2)
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")


#Set settings for Heston
settings <- sim.setup(mat=6.5/(24*7*52), Npath = 2, omega = 1.6*10^(-5)) #6.5 hours

#Get results with and without bursts. OBS: I have replaced c_1=3 with c_1 = 0.3
Heston <- sim.heston(settings)
Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = 1.5, beta = 0.1,reverse = F)
Heston_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05,c_1 = 0.016,alpha=0.8, reverse = F)
Heston_jump <- sim.addjump(Heston,burst_time=0.5,interval_length=0.05, c_1 = 0.016, alpha = 0.8)


#Get a single path
path <- 2

Heston_path <- sim.path(path,Heston)$Y
vb_path <- sim.path(path,Heston_vb)$Y
vbdb_path <- sim.path(path,Heston_vbdb)$Y
jump_path <- sim.path(path, Heston_jump)$Y

#Set time index for plot
time_index <- Heston$time/settings$mat


#Only plot from 0.35 to 0.65
index <- ((time_index >=0.35) & (time_index <=0.65))
Heston_path <- Heston_path[index]
vb_path <- vb_path[index]
vbdb_path <- vbdb_path[index]
jump_path <- jump_path[index]
time_index <- time_index[index]


#----------------------------- NICE PLOT (NIELS' REGRESSION CODE) -----------------------------  

df_y_values <- data.frame(cbind(Heston_path,vb_path, vbdb_path, jump_path))
names <- c("No burst","+Volatility burst", "+drift burst", "+jump")

tmp <- list()
for (i in 1:4) {
  tmp[[i]] <- data.frame(
    x_akse = time_index,
    y_akse = df_y_values[,i],
    rx = names[i]
  )
}
tmp <- do.call(rbind, tmp)

ggplot(tmp, aes(x_akse, y_akse,color = rx)) +
  geom_line() +
  scale_color_manual("", values = c("black", "red", "blue", "purple")) +
  xlab("time") + ylab("log-return")  + 
  #theme(legend.position = c(1,0),
  #      legend.justification = c(1,0),
  #      plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Log-return of asset") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))
