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
hd <- 2*300 #(seconds)
hv <- hd                                    # ATTEMPT RATIO of 1 FOR DB
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
SPY <- data.TtoStar(SPY.T, "day", 0.95)
bitfinex<-data.TtoStar(bitfinex.T, "day", 0.95)
bitmex<-data.TtoStar(bitmex.T, "day", 0.95)
kraken<-data.TtoStar(kraken.T, "day", 0.95)

# TABLE
if(F){
  length(unique(SPY$day))
  sum(SPY$db_0.95)
  sum(SPY$Tstar>5)
  sum(SPY$Tstar>6)
  sum(SPY$Tstar>7)
  sum(SPY$Tstar>10)
  
  length(unique(bitfinex$day))
  sum(bitfinex$db_0.95)
  sum(bitfinex$Tstar>5)
  sum(bitfinex$Tstar>6)
  sum(bitfinex$Tstar>7)
  sum(bitfinex$Tstar>10)
  
  length(unique(bitmex$day))
  sum(bitmex$db_0.95)
  sum(bitmex$Tstar>5)
  sum(bitmex$Tstar>6)
  sum(bitmex$Tstar>7)
  sum(bitmex$Tstar>10)
  
  length(unique(kraken$day))
  sum(kraken$db_0.95)
  sum(kraken$Tstar>5)
  sum(kraken$Tstar>6)
  sum(kraken$Tstar>7)
  sum(kraken$Tstar>10)
}

# PLOTTERIA
require(ggplot2)
df <- data.frame(type = factor(c(rep("SPY",length(SPY$Tstar)), rep("Bitfinex",length(bitfinex$Tstar)),
                                 rep("BitMEX",length(bitmex$Tstar)), rep("Kraken",length(kraken$Tstar)))),
                 Tstar = c(SPY$Tstar, bitfinex$Tstar, bitmex$Tstar, kraken$Tstar))

ggplot(df, aes(x = Tstar, fill = type)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 4) +
  xlab("Value of T* (daily)") +
  ylab("Density")


