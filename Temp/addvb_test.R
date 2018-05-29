set.seed(100)
library(microbenchmark)
setwd(Sys.getenv("masters-thesis"))
source("Simulation/Heston.R")
source("Simulation/Bursts.R")

#################### PARAMETERS ####################
omega2 <- 2.64*10^(-10)
omega <- sqrt(omega2)
K2 <- 0.5 #K2
n <- 23400
mat <- 6.5/(24*7*52)
dt <- mat/n
Npaths <- 1000 #One year of S&P-500
sigma2 <- 0.0457
sigma <- sqrt(sigma2)
(noise_ratio <- omega/(sqrt(dt)*sigma)) #10.66

### SIMULATE
settings <- sim.setup(mat=mat, Npath = Npaths, Nsteps = n, omega = omega) #6.5 hours
heston <- sim.heston(settings)

source("Simulation/Bursts.R")
p0 <- Sys.time()
test1 <- sim.addvb(heston,c_2 = 5.44*10^(-4), beta = 0.45)
print(Sys.time()-p0)

p0 <- Sys.time()
test2 <- sim.addvb.2.0(heston,c_2 = 5.44*10^(-4), beta = 0.45)
print(Sys.time()-p0)

#No difference
max(abs(test1$X-test2$X))
max(abs(test1$Y-test2$Y))
max(abs(test1$vol-test2$vol))
max(abs(test1$time-test2$time))

#Test speed up
times <- microbenchmark(sim.addvb(heston,c_2 = 5.44*10^(-4), beta = 0.45),
                        sim.addvb.2.0(heston,c_2 = 5.44*10^(-4), beta = 0.45),
                        times = 1)

paste("Relative speed-up:",round(max(times$time)/min(times$time),2))
