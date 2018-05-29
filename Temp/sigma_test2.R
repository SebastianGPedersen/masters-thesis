set.seed(100)
library(ggplot2)
library(latex2exp)
library(microbenchmark)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")
source("Simulation/Jumps.R")
source("Estimation/pre-average.R")
source("Estimation/estimates.R")
source("Estimation/estimates_reloaded.R")
source("Estimation/estimates_revolution.R")
source("Estimation/gumbel.R")
source("Kernels/kernels.R")

#################### PARAMETERS ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 100 #One year of S&P-500
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

sigma_2 <- est.sigma.mat.2.0(data = path, hv = h_mu, lag = lag)$sig[,(lag+1):n]
sigma_3 <- est.sigma.mat.3.0(data = path, hv = h_mu, lag = lag)$sig[,(lag+1):n]

times <- microbenchmark(est.sigma.mat.2.0(data = path, hv = h_mu, lag = lag),
                        est.sigma.mat.3.0(data = path, hv = h_mu, lag = lag),
                        times = 1)

paste("Maximum difference across all estimators:",max(abs(sigma_2/sigma_3-1))) #Max difference between the two
paste("Relative speed-up:",round(max(times$time)/min(times$time),0))
