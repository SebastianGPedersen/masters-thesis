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

require(ggplot2)
# REAL DATA
sp <- data.getFull()

sp<-data.dayID(sp)

start <- exp(sp$logPrice[1])

close <-exp(sp$logPrice[length(sp$logPrice)])

return <- sp$logPrice[length(sp$logPrice)]-sp$logPrice[1]

high <- max(exp(sp$logPrice))

low <- min(exp(sp$logPrice))

#only every 1k-ths observation

plot.data <- data.table(DateTime = sp$DateTime[seq(1,length(sp$DateTime), by = 1000)], Price = exp(sp$logPrice[seq(1,length(sp$DateTime), by = 1000)]))

g1 <- ggplot(data = plot.data) +
  geom_line(aes(x = DateTime, y = Price, group = 1, col = "SPY")) +
  ylab("Price") +
  xlab("Time")
  #scale_x_datetime("Time")
g1

rm(sp)
# BIT PENGE

bitfinex <- data.getbitdata("bitfinex")
bitmex <- data.getbitdata("bitmex")
kraken <- data.getbitdata("kraken")

# FINEX
(start <- exp(bitfinex$logPrice[1]))
(close <- exp(bitfinex$logPrice[length(bitfinex$logPrice)]))

(return <- bitfinex$logPrice[length(bitfinex$logPrice)]-bitfinex$logPrice[1])

(high <- max(exp(bitfinex$logPrice)))

(low <- min(exp(bitfinex$logPrice)))

# MEX
(start <- exp(bitmex$logPrice[1]))
(close <- exp(bitmex$logPrice[length(bitmex$logPrice)]))

(return <- bitmex$logPrice[length(bitmex$logPrice)]-bitmex$logPrice[1])

(high <- max(exp(bitmex$logPrice)))

(low <- min(exp(bitmex$logPrice)))

# Kraken
(start <- exp(kraken$logPrice[1]))
(close <- exp(kraken$logPrice[length(kraken$logPrice)-1]))

(return <- kraken$logPrice[length(kraken$logPrice)-1]-kraken$logPrice[1])

(high <- max(exp(kraken$logPrice)))

(low <- min(exp(kraken$logPrice)))

# PLOT

plot.data1 <- data.frame(DateTime_fin = bitfinex$DateTime[seq(1,length(bitfinex$DateTime), by = 1000)], Price_fin = exp(bitfinex$logPrice[seq(1,length(bitfinex$logPrice), by = 1000)]))
plot.data2 <- data.frame(DateTime_mex = bitmex$DateTime[seq(1,length(bitmex$DateTime), by = 1000)], Price_mex = exp(bitmex$logPrice[seq(1,length(bitmex$logPrice), by = 1000)]))
plot.data3 <- data.frame(DateTime_kraken = kraken$DateTime[seq(1,length(kraken$DateTime), by = 1000)], Price_kraken = exp(kraken$logPrice[seq(1,length(kraken$logPrice), by = 1000)]))

g1 <- ggplot() +
  geom_line(data = plot.data1, aes(x = DateTime_fin, y = Price_fin, group = 1, col = "BitFinex")) +
  geom_line(data = plot.data2, aes(x = DateTime_mex, y = Price_mex, group = 1, col = "BitMEX")) +
  geom_line(data = plot.data3, aes(x = DateTime_kraken, y = Price_kraken, group = 1, col = "Kraken")) +
  ylab("Price") +
  xlab("Time")
#scale_x_datetime("Time")
g1

plot(tail(kraken$DateTime), type = "l")

bitfinex$Date[1]
bitfinex$Date[length(bitfinex$Date)]

bitmex$Date[1]
bitmex$Date[length(bitmex$Date)]

kraken$Date[1]
kraken$Date[length(kraken$Date)]
