DT <- data.getbitdata("bitfinex")
DT <- data.getbitdata("kraken")
DT <- data.getbitdata("bitMEX")

# plot(diff(DT$DateTime))
# 
# dt <- diff(DT$Time)
# plot((which(dt>=60)))
# density(dt[(which(dt>=60))])
# 
# ddt <- diff((which(dt>=60)))
# which(ddt == max(ddt))
# ddt[which(ddt == max(ddt))-1]


t <- DT$Time
# t <-c( 1, 2, 100 , 101 , 103, 104, 105, 106, 107, 1000, 1001, 1002, 1003, 1500, 1501)
maxTol <- 120
dt <- diff(t)

# t[ which(dt>=maxTol)]

# which(dt>=maxTol)
ddt <- diff(c(0,which(dt>=maxTol), length(t)))
# ddt

# c(0,which(dt>=maxTol), length(t))
# ddt <- c(1,2,3)
crit <- ddt==max(ddt[ddt!=max(ddt)]) # second longest interval
crit <- ddt==max(ddt) # longest interval. 
maxPeriod <- which(crit)

StartEnd <-  c(c(0,which(dt>=maxTol), length(t))[maxPeriod] + 1 , c(0,which(dt>=maxTol), length(t))[maxPeriod+1])
DT2 <- DT[StartEnd[1]:StartEnd[2],]

DT2

quantile(diff(DT2$Time), probs = c(75:100/100))
quantile(diff(DT2$Time), probs = c(990:1000/1000))
quantile(diff(DT2$Time), probs = c(9990:10000/10000))
quantile(diff(DT2$Time), probs = c(99990:100000/100000))
quantile(diff(DT2$Time), probs = c(999990:1000000/1000000))
quantile(diff(DT2$Time), probs = c(9999990:10000000/10000000))
plot(9999000:10000000/10000000, quantile(diff(DT2$Time), probs = c(9999000:10000000/10000000)))

