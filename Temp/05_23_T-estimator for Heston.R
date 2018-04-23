library(ggplot2)
source(paste(Sys.getenv("masters-thesis"),"Simulation/Heston.R",sep="/"))
source(paste(Sys.getenv("masters-thesis"),"Simulation/Bursts.R",sep="/"))

#################### PARAMETERS THAT DON'T CHANGE ####################
omega <- 0.0000225
omega2 <- omega^2
ksq <- 0.5 # K2
phi_1 <- 1 #int(g'(x)^2)
phi_2 <- 1/12 #int(g^2)
int_g <- 1/4 #int(g)


#################### PARAMETERS CHANGING WITH N ####################
n <- 1000
k_n <- floor(1/2*n^(1/2))
hd <- floor(10*n^(-1/4))
  
  
#################### Simulation ####################
Npath <- 500
settings <- sim.setup(mat=6.5/(24*7*52), Npath = Npath, Nsteps = n, omega = 0.0000225) #6.5 hours
Heston <- sim.heston(settings)


#################### Pre-averaging ####################

#We need to transpose Y for dy to work properly
Y <- t(as.matrix(Heston$Y))
dy <- diff(Y)

#Pre-avg (can't take vector)
k_n <- floor(1/2*sqrt(length(Heston$time)))
pre_y <- matrix(NA, nrow = Npath, ncol = n-k_n)

for (i in 1:Npath){
  pre_y[i,] <- est.NewPreAverage(dy[,i],k_n)
}
  
#################### T estimator ####################

#Calculate T estimator
mu_hat <- 1/(h_n*k_n*int_g)*sum(pre_y)
sigma_hat <- 
T_hat <- 

#################### PLOT ####################

#Et loop mere senere til at vise konvergens

Heston_vb <- sim.addvb(Heston,burst_time = 0.5, interval_length = 0.05, c_2 = 0.1, beta = 0.1)
Heston_vbdb <- sim.adddb(Heston_vb, burst_time=0.5,interval_length=0.05,c_1 = 0.03,alpha=0.8)

#Get a single path
path <- 2

Heston_path <- sim.path(path,Heston)$Y
vb_path <- sim.path(path,Heston_vb)$Y
vbdb_path <- sim.path(path,Heston_vbdb)$Y

#Set time index for plot
time_index <- Heston$time/settings$mat


#Only plot from 0.35 to 0.65
index <- ((time_index >=0.35) & (time_index <=0.65))
Heston_path <- Heston_path[index]
vb_path <- vb_path[index]
vbdb_path <- vbdb_path[index]
time_index <- time_index[index]


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
