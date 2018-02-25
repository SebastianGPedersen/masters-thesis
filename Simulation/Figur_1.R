#Module takes 13 seconds to run

library(ggplot2)
source(paste(Sys.getenv("masters-thesis"),"Simulation/Heston.R",sep="/"))
source(paste(Sys.getenv("masters-thesis"),"Simulation/Bursts.R",sep="/"))


#Set settings for Heston
settings <- sim.setup(mat=1)

#Get results with and without bursts
Heston <- sim.heston(settings)
Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = 0.15, beta = 0.4)
Heston_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05,c_1=3,alpha=0.75)



#Get a single path
path = 2

Heston_path = sim.path(path,Heston)$Y
vb_path = sim.path(path,Heston_vb)$Y
vbdb_path = sim.path(path,Heston_vbdb)$Y

#Set time index for plot
time_index = Heston$time/settings$mat


#Only plot from 0.35 to 0.65
index = ((time_index >=0.35) & (time_index <=0.65))
Heston_path = Heston_path[index]
vb_path = vb_path[index]
vbdb_path = vbdb_path[index]
time_index = time_index[index]


#----------------------------- NICE PLOT (NIELS' REGRESSION CODE) -----------------------------  

df_y_values <- data.frame(cbind(Heston_path,vb_path, vbdb_path))
names <- c("No burst","+Volatility burst", "+drift burst")

tmp <- list()
for (i in 1:3) {
  tmp[[i]] <- data.frame(
    x_akse = time_index,
    y_akse = df_y_values[,i],
    rx = names[i]
  )
}
tmp <- do.call(rbind, tmp)

ggplot(tmp, aes(x_akse, y_akse,color = rx)) +
  geom_line() +
  scale_color_manual("", values = c("black", "red", "blue")) +
  xlab("time") + ylab("log-return")  + 
  #theme(legend.position = c(1,0),
  #      legend.justification = c(1,0),
  #      plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Log-return of asset") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

