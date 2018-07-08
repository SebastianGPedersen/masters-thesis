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

data <-data.getFull()

# MARK DAYS
data <- data.dayID(data)

# ESTIMATION PARAMETERS
hd <- 300 #(seconds)
hv <- 12*hd
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)

#
Tdata <- data.TforId(data, "day", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)

Tstar_daily <- data.TtoStar(Tdata, "day", 0.95)

# choose bursts = T
bursts <- Tstar_daily[db_0.95 == T,]

# PICK ONE OUT
(burst <- bursts[15,])

# PLOTTERIA
data.plot_db(data, burst$DateTime, hd = hd, hv = hv, lag = 10)

weekdays(as.Date(burst$DateTime))

# PLOT ZOOMED OUT

fullday <- data[day == burst$day,]


require(ggplot2)
plotdata <- data.frame(Date = fullday$DateTime, Price = exp(fullday$logPrice))

ggplot(data = plotdata, aes(x = Date)) +
    geom_line(aes(y = Price, colour = "Price"))

fullday <- data[day == burst$day,]

poi <- as.POSIXct("2014-03-19 15:00:00", tz = "UTC")

intra <- data[DateTime > poi-20*60 & DateTime < poi+20*60, ]

plotdata <- data.frame(Date = intra$DateTime, Price = exp(intra$logPrice))

ggplot(data = plotdata, aes(x = Date)) +
  geom_line(aes(y = Price, colour = "Price"))
