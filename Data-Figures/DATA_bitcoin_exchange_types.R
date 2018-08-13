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
hv <- 12*hd                                          # TRY RATIO 1
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)

# GET DATA AND ROLL
bitfinex <- data.getbitdata("bitfinex")
bitmex <- data.getbitdata("bitmex")
kraken <- data.getbitdata("kraken")

poi <- as.POSIXct("2018-06-10 12:00:00", tz = "UTC")

window = 5*60

finex <- bitfinex[DateTime < poi + window & DateTime > poi - window, ]
mex <- bitmex[DateTime < poi + window & DateTime > poi - window, ]
krak <- kraken[DateTime < poi + window & DateTime > poi - window, ]

# PLOTTERIA

ggplot() +
  geom_line(aes(x = DateTime, y = exp(logPrice), colour = "BitFinex"), data = finex, size = 1) + 
  geom_line(aes(x = DateTime, y = exp(logPrice), colour = "BitMEX"), data = mex, size = 1) + 
  geom_line(aes(x = DateTime, y = exp(logPrice), colour = "Kraken"), data = krak, size = 1) + 
  xlab("Time") + ylab("Price")
  
