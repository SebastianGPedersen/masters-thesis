setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
source("spy/dataFunctions.R")
source("spy/datahandling.R")

# APPLY ON DAY #

data<-data.getFull()
data<-data.dayID(data, id = "day")


dt <- data[day == 1, ] #tmp
burst_time <- data.TstarforId(data = dt, id_name = "day", id_number = 1, hd = 1, t.freq = 10, lag = 10, offset = 10, conf = 95)$DateTime
window <- 20*60 # 20-minutes


data.plot_db<-function(data, burst_time, window = 20*60, hd, hv, lag = 10){
  # NEEDS SCALE/OFFSET PARAMETERS FOR AESTETICS
  data <- dt
  # FIND WINDOW
  start<- burst_time-window
  end <- burst_time+window
  start_index <- data.date_To_tindex(start, dt$DateTime)
  end_index <- data.date_To_tindex(end, dt$DateTime)
  
  # EXTRACT DATA
  focus <- dt[start_index:end_index, ]
  
  price <- exp(focus$logPrice)
  Y <- diff(dt$logPrice) # we need all data for estimation
  time <- diff(dt$Time)
  DateTime <- focus$DateTime
  
  # ESTIMATE T
  est_data <- list(time = time, Y = Y)
  t.freq <- start_index:end_index       # check this
  Ts <- arima.sim(model=list(ar=0.9), n = length(start_index:end_index))/2 #$test      # change this to actual T
  
  # PLOTTERIA
  require(ggplot2)
  plotdata <- data.frame(Date = DateTime, Price = price, T.statistic = Ts) # mange punkter uden ændring
  
  offset <- mean(price)
  scale <- (max(Ts)-min(Ts))/(max(price)-min(price)) # check this later
  p <- ggplot(data = plotdata, aes(x = Date)) +
    geom_area(aes(y = T.statistic)) + #rescale this
    geom_line(aes(y = (Price-offset)*scale, colour = "Price"), size = 1) + # size might change
    scale_y_continuous(sec.axis = sec_axis(trans = ~./scale+offset, name = "Price"))
  p
  return(p)
}

data.plot_db(data[day == 1, ], burst_time = burst_time, window = window, hd = 1, hv = 1, lag = 10)
