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
hv <- 12*hd                                               
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)

# GET DATA AND ROLL
SPY <- data.getFull()
bitfinex <- data.getbitdata("bitfinex")
bitmex <- data.getbitdata("bitmex")
kraken <- data.getbitdata("kraken")

bitfinex <- bitfinex[, all := 1]
bitmex <- bitmex[, all := 1]
kraken <- kraken[, all := 1]

bitfinex.T<-data.TforId(bitfinex, "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
bitmex.T<-data.TforId(bitmex, "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
kraken.T<-data.TforId(kraken, "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)

# MARK DAYS
SPY <- data.dayID(SPY)
bitfinex.T <- data.dayID(bitfinex.T)
bitmex.T <- data.dayID(bitmex.T)
kraken.T <- data.dayID(kraken.T)

SPY.T <- data.TforId(SPY, "day", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)

gc()

# TSTARIA
SPY.Tstar <- data.TtoStar(SPY.T, "day", 0.95)
bitfinex.Tstar<-data.TtoStar(bitfinex.T, "day", 0.95)
bitmex.Tstar<-data.TtoStar(bitmex.T, "day", 0.95)
kraken.Tstar<-data.TtoStar(kraken.T, "day", 0.95)

# ##########################################
#                                          #
#               PLOTS - BITFINEX           #
#                                          #
# ##########################################

bitfinex.bursts <- bitfinex.Tstar[db_0.95 == T,]

# PICK ONE OUT
bitfinex.burst <- bitfinex.bursts[77,]

# PLOTTERIA
data.plot_db(bitfinex, bitfinex.burst$DateTime, hd = hd, hv = hv, lag = 10, blue = T)

# FIND WORST
min(bitfinex.T$Tval)

# ##########################################
#                                          #
#               PLOTS - BITMEX             #
#                                          #
# ##########################################

bitmex.bursts <- bitmex.Tstar[db_0.95 == T,]

# PICK ONE OUT
bitmex.burst <- bitmex.bursts[77,]

# PLOTTERIA
data.plot_db(bitmex, bitmex.burst$DateTime, hd = hd, hv = hv, lag = 10, blue =T)

# ##########################################
#                                          #
#               PLOTS - KRAKEN             #
#                                          #
# ##########################################

kraken.bursts <- kraken.Tstar[db_0.95 == T,]

# PICK ONE OUT
kraken.burst <- kraken.bursts[62,]

# PLOTTERIA
data.plot_db(kraken, kraken.burst$DateTime, hd = hd, hv = hv, lag = 10, blue =T)
