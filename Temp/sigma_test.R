set.seed(100)
library(ggplot2)
library(latex2exp)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")
source("Estimation/estimates_reloaded.R")
source("Estimation/gumbel.R")
p0 <- Sys.time()

#################### PARAMETERS ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 10 #One year of S&P-500
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
lag <- 10 #Temporary. Should be 100
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66
h_mu <- 5 / (60*24*7*52) #5min bandwidth of drift-estimator
T_interval <- 5 / (60*60*24*7*52) #5sec between T_calculations
#T_interval <- 1 / (60*24*7*52) #1min between T_calculations
m <- mat/T_interval
threshold <- q95(m)
#threshold <- qnorm(0.975)


### SIMULATE
settings <- sim.setup(mat=mat, Npath = Npaths, Nsteps = n, omega = omega) #6.5 hours
path <- sim.heston(settings)
#Transform Y to dY in path$Y
path$Y <- t(diff(t(as.matrix(path$Y))))
#Extract the relevant part
path$X <- NULL
path$vol <- NULL

### Indexes from sigma-estimator
desired_indices <- floor(T_interval/dt)*(1:(m-1))[-(1:2)]

### Compare calculations
path_single <- path
path_single$Y <- path_single$Y[1,]


p0 <- Sys.time()
sigma_1 <- est.sigma.mat.next(data = path, hv = h_mu, t.index = desired_indices, lag = lag)$sig
time1 <- as.numeric(difftime(Sys.time(),p0,units = "secs"))

p0 <- Sys.time()
sigma_2 <- est.sigma.mat.2.0(data = path, hv = h_mu, lag = lag)$sig[,desired_indices]
(time2 <- as.numeric(difftime(Sys.time(),p0,units = "secs")))

paste("Maximum difference across all estimators:",max(abs(sigma_2/sigma_1-1))) #Max difference between the two
paste("Relative speed-up:",round(as.numeric(time1/time2),0))

plot(sigma_1[1,] - sigma_2[1,], type="l")
