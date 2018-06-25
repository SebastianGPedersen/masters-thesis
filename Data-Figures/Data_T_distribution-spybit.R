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
hd <- 600 #(seconds)
hv <- 12*hd                                             # ATTEMPT WITH BANDWIDTH RATIO 1
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

gc()

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

# DISTRIBUTION OF T
SPY.T<-SPY.T$Tval
bitfinex.T<-bitfinex.T$Tval
bitmex.T<-bitmex.T$Tval
kraken.T<-kraken.T$Tval

gc()

# PREP TABLE
require(moments)
SPY<-c(mean(SPY.T), sd(SPY.T), kurtosis(SPY.T))
names(SPY) <- c("mean", "sd", "kurtosis")
print(SPY)

bitfinex<-c(mean(bitfinex.T), sd(bitfinex.T), kurtosis(bitfinex.T))
names(bitfinex) <- c("mean", "sd", "kurtosis")
print(bitfinex)

bitmex<-c(mean(bitmex.T), sd(bitmex.T), kurtosis(bitmex.T))
names(bitmex) <- c("mean", "sd", "kurtosis")
print(bitmex)

kraken<-c(mean(kraken.T), sd(kraken.T), kurtosis(kraken.T))
names(kraken) <- c("mean", "sd", "kurtosis")
print(kraken)



# PLOT DENSITY
require(ggplot2)
df <- data.frame(type = factor(c(rep("SPY",length(SPY.T)), rep("Bitfinex",length(bitfinex.T)),
                                 rep("BitMEX",length(bitmex.T)), rep("Kraken",length(kraken.T))
)),
Tvalue = c(SPY.T, bitfinex.T, bitmex.T, kraken.T))

ggplot(df, aes(x = Tvalue, fill = type)) +
  geom_density(alpha = 0.5) +
  xlab("Value of T-estimator") +
  ylab("Density")



