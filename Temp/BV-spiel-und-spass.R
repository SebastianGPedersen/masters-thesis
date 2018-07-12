setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/BlackScholes.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("kernels/kernels.R")
source("vol/vol_estimators.R")

# FACT OR FRICTION TABLE 2 #

# EXPERIMENT WITH INTEGRATED VOLATILITY AND BV #

# K = theta*sqrt(N)+theta*sqrt(N)%%2 # FORCE EVEN

# ~~~BS~~~
drift <- 0.0; vol <- 0.2; noise <- 0#1.6*10^-5
mat <- 6.5/(52*7*24)
Nsteps <- 50000; Npath <- 3

data <- sim.BlackScholes(drift, vol, noise, mat, Nsteps, Npath)

x <- data$Y[1,]

R <- diff(x)

# Realized approx
sum(R^2)
# Actual
vol^2*mat
# BV
vol.est.BV(R)

# ~~~HESTON~~~
hest<-sim.heston(sim.setup(Npath = Npath, Nsteps = Nsteps))

x<-hest$X[1,]

R <- diff(x)

# Realized approx
sum(R^2)
# Actual
sum(hest$vol[1,]*diff(hest$time)[1])
# BV
vol.est.BV(R)

# VB
burst_time <- 0.5; interval_length <- 0.05; alpha <- 0.8; beta <- 0.1;
c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha); c_2 <- sqrt((1-2*beta)*0.001^2/(10/(60*24*7*52))^(1-2*beta));
length(hest)
hestvb <- sim.addvb(Heston_res = hest, burst_time = burst_time, interval_length = interval_length,
                    c_2 = c_2, beta = beta, reverse = F, recenter = F)

plot(hestvb$vol[1,], type = "l")

x<-hestvb$X[1,]

R <- diff(x)

# Realized approx
sum(R^2)
# Actual
sum(hest$vol[1,]*diff(hest$time)[1])
# BV
vol.est.BV(R)

