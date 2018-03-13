setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("kernels/kernels.R")

bursts <- function(settings, burstsetting){
  alpha <- burstsetting$alpha
  beta <- burstsetting$beta
  burst_time <- burstsetting$burst_time
  interval_length <- burstsetting$interval_length
  c_1 <- burstsetting$c_1
  c_2 <- burstsetting$c_2
  
  sims<-sim.heston(setting)
  
  sims.vb<-sim.addvb(sims,    burst_time = burst_time, interval_length = interval_length,
                                c_2 = c_2, beta  = beta)
  
  sims.db<-sim.adddb(sims.vb, burst_time = burst_time, interval_length = interval_length,
                                c_1 = c_1,  alpha = alpha)
                    

  #plot(sim.db$Y, type = "l", col = "blue")
  #lines(sim.vb$Y, col = "red")
  #lines(sim$Y)
  
  #Get a single path
  path = 1
  
  Heston_path = sim.path(path,sims)$Y
  vb_path = sim.path(path, sims.vb)$Y
  vbdb_path = sim.path(path, sims.db)$Y
  
  #Set time index for plot
  time_index = sims$time/settings$mat
  
  
  #Only plot from 0.35 to 0.65
  index = ((time_index >=0.35) & (time_index <=0.65))
  Heston_path = Heston_path[index]
  vb_path = vb_path[index]
  vbdb_path = vbdb_path[index]
  time_index = time_index[index]
  
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
    ggtitle("Log-return of asset") +
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  #return(sim$Y[11700]-sim.db$Y[11700])
}

require(ggplot2)
setting <- sim.setup(Npath = 2, Nstep = 23400, omega = 0.0225/10000)
burst<-sim.burstsetting(alpha = 0.55, beta = 0.2 ,c_1 = 0.2, c_2 = 0.03)

set.seed(1234)

bursts(setting, burst)
