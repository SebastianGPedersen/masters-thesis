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
firstDay<-selectDays(fullData, as.Date("2014-06-05"), nDays = 2)

# A QUICK VIEW OF DATA
plot(firstDay$logPrice, type="l")

# GLUE DATA TO FRAME AND FIND TIME INDEX FOR ESTIMATION
data <- list(time = firstDay$Time, Y = firstDay$logPrice)
tind<-timePoints(firstDay, timeOffset = 120, initialDelay = 600)
tind<-tind[1:(length(tind)-1)] # We dont know DY in the latest, do we?

# PRE AVERAGE
horizon<-floor(1*sqrt(length(data$time))) # FORCE INTEGER
k = horizon + (horizon%%2) # FORCE EVEN
prev <- c(rep(0,k-2),est.PreAverage(data$Y, k)) # Kim does not divide by k - scaling does not change results
data <- list(time = data$time, Y = prev)

# BANDWIDTH
dt <- diff(data$time)[1]
hd <- 300*dt
hv <- 1500*dt

conf = 0.95

# Estimation of mu/sig
mu<-est.mu.next(data = data, hd = hd, t.index = tind)

sig <- est.sigma.next(data, hv=hv, t.index = tind, lag = 15)#"auto")

# Calculate T
Tstat<-teststat(mu, sig, hd, hv)

# Calculate T*
Tstar<-tstar(Tstat)$tstar

# fit rho
rho <- est.rho(Tstat$test)

z<-est.z_quantile(rho$m, rho$rho, conf)$qZm
res<-Tstar>=z

a<-c(res, Tstar, z)
names(a) <- c("T/F", "Tstar", "z")

print(a)

plot(Tstat$test, type = "l")

mean(Tstat$test)
var(Tstat$test)
qqnorm(Tstat$test)
qqline(Tstat$test)
