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
source("vol/vol_estimators.R")

# ESTIMATION PARAMETERS
hd <- 300 #(seconds)
hv <- 12*hd                                          # TRY RATIO 1
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)

# GET DATA AND ROLL
bitfinex <- data.getbitdata("bitfinex")

bitfinex <- bitfinex[, all := 1]

bitfinex.T<-data.TforId(bitfinex, "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
# MARK DAYS
bitfinex.T <- data.dayID(bitfinex.T)

gc()

# TSTARIA
bitfinex.Tstar<-data.TtoStar(bitfinex.T, "day", 0.95)

BV<-vol.est.IVestPathwise()

