setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/rescaling.R")
source("kernels/kernels.R")
source("spy/dataFunctions.R")
source("spy/datahandling.R")

# ESTIMATION PARAMETERS
hd <- 300 #(seconds)
hv <- hd*12                                   # ATTEMPT RATIO of 1 FOR DB
lag = 30
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)

# GET DATA AND ROLL
SPY <- data.getFull()

# MARK DAYS
SPY <- data.dayID(SPY)

SPY.T <- data.TforId(SPY, "day", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)

gc()

# TSTARIA
SPY <- data.TtoStar(SPY.T, "day", 0.95)
burst <- SPY[db_0.95 == T, ]


t <- strftime(burst$DateTime, format = "%H:%M:%S", tz = "UTC")

# Set all to same date
ntimes.unsorted <- as.POSIXct(t, format = "%H:%M:%S", tz = "UTC")

# SORT
pre.sort <-data.table(Numeric = as.numeric(ntimes.unsorted), DateTime = ntimes.unsorted, absT = burst$Tstar)
sorted <- setkey(pre.sort)

ntimes<-sorted$DateTime

data <- data.frame(times = Times)
start <- as.POSIXct("2018-08-02 09:30:00 UTC")
end <-  as.POSIXct("2018-08-02 16:00:00 UTC")
breaks <- seq(from = start, to = end, by =  60*60)

# LAV DET MANUELT

require(ggplot2)

ggplot(data = data) +
  geom_histogram(aes(ntimes)) +
  xlab("Time") + ylab("Density")

