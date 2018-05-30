set.seed(100)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Estimation/estimates.R")
source("Estimation/estimates_reloaded.R")

#################### PARAMETERS ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 2 #One year of S&P-500
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66
h_mu <- 5 / (60*24*7*52) #5min bandwidth of drift-estimator

### SIMULATE
settings <- sim.setup(mat=mat, Npath = Npaths, Nsteps = n, omega = omega) #6.5 hours
path <- sim.heston(settings)

#Transform Y to dY in path$Y
path$Y <- t(diff(t(as.matrix(path$Y))))

#Extract the relevant part
path$X <- NULL
path$vol <- NULL

### Compare calculations
path_single <- path
path_single$Y <- path_single$Y[1,]


mu <- est.mu(data=path_single,t.index = 1, hd = h_mu)$mu
mu2 <- est.mu.mat.2.0(data=path,hd = h_mu)$mu[1,1]

#correct mu
dy <- path$Y[1,1]
kernel <- kern.leftexp$kern((path$time[1]-path$time[2])/h_mu)
mu_hat <- 1/h_mu*kernel*dy
mu_hat_wrong_index <- 1/h_mu*dy
