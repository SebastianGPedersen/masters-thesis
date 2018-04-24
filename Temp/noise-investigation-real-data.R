setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/pre-average.R")
source("kernels/kernels.R")
source("SPY/datafunctions.R")

# IMPORT REAL DATA
fullData<-readRDS(paste0(Sys.getenv("masters-thesis-data"),"/SPY/2014_SPY_Vol_Avg.rds"))

firstDay<-selectDays(fullData, as.Date("2014-06-05"), nDays = 1)

# A QUICK VIEW OF DATA
plot(firstDay$logPrice, type="l")

#qqnorm(diff(firstDay$logPrice))

# GLUE DATA TO FRAME AND FIND TIME INDEX FOR ESTIMATION
data <- list(time = firstDay$Time, Y = firstDay$logPrice)
tind<-timePoints(firstDay, timeOffset = 120, initialDelay = 600)
tind<-tind[1:(length(tind)-1)] # We dont know DY in the latest

theta = 1
# PRE AVERAGE
if(1 == 1){
  horizon<-floor(theta*sqrt(length(data$time))) # FORCE INTEGER
  k = horizon + (horizon%%2) # FORCE EVEN
  prev <- c(rep(0,k-2),est.PreAverage(data$Y, k)) # Kim does not divide by k - scaling does not change results
}

# No preaverage
if(0==1){
  k<-1
  prev <- data$Y; raw = data$Y
}

data <- list(time = data$time, Y = prev, raw = data$Y)
# BANDWIDTH
dt <- diff(data$time)[1]
hd <- 1500*dt
hv <- 1500*dt

# Estimation of noise
mu <- est.mu.new(data = data, hd = hd, t.index = tind, kn = 1)
sig <- est.sigma.new(data = data, hv = hv, t.index = tind, noisefun = est.noise.iid.next,
                       kn = k, theta = theta)
mean(sig$noise)
plot(sig$noise, type = "l")


