setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")
source("Estimation/estimates.R")

# Time parameters
mat <- 1
t <- 0:10000  # time_points
dt <- mat/t[length(t)]

# Define sigma
sig <- 0.1
sig2 <- sig^2
omega2 <- 1
omega <- 1
ksq <- 0.5 # K2
hd <- 0.01 #bandwidth in mu
hv <- 0.05 #bandwidth in sigma


# Simulate a Brownian motion (X) and one with mcrostrucutre noise (Y)
x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
eps <- rnorm(length(t), mean = 0, sd = sqrt(omega))
x <- cumsum(x)
y <- x + eps


# Plot part of x and y to validate it is done correctly
plot(x[1:100],type="l")
plot(y[1:100],type = "l")
#Støjen overgår tydeligt X som den burde



############ Test convergence of mu without noise (eq. 12 and 13) ############

N <- 5000
mu <- numeric(N)

for(i in 1:N) #Takes approx 10sec
{
  #simulate brownian motion
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x)
  
  #put in data.frame for use in est.mu function
  dataX<-data.frame(time = t*dt, Y = x)
  
  mu[i] <- est.mu(dataX, hd, kern.leftexp, t.index = 500)$mu[1] #Take out a random timeindex (500)
}

mu <- sqrt(hd)*mu #This should be N(0,ksq*omega^2) distributed from eq. (13)

#Normalize
mu_normal <- mu*sqrt(1/(ksq*sig2))

#test - it all looks really fine
mean(mu_normal)
var(mu_normal)
plot(mu_normal)
qqnorm(mu_normal)




############ Test convergence of sigma without noise (convergence of eq. 14 to sigma) ############

#Simulate the BM
x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
x <- cumsum(x)

#Put in dataframe
dataX<-data.frame(time = t*dt, Y = x)

#Get sigma estimates on 1000 different time_points (vigtigt at bandwidth er rigtig lille)
sigma_1000_times <-est.sigma.raw(data = dataX, hv, kern.leftexp, t.index = seq(1000, 10000, by = 10))$sigRaw

est_sig2 <- sigma_1000_times^2

#Test that this is close to sig2 (0.01)
plot(est_sig2)

#It looks centralized plus minus 20% so not bad




############ Test convergence of mu with noise (convergence of eq. 25 to Theorem 5 p. 43) ############

N <- 5000
mu <- numeric(N)

for(i in 1:N) #Takes approx 30sec
{
  #simulate brownian motion with noise
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x)
  eps <- rnorm(length(t), mean = 0, sd = sqrt(omega2))
  y <- x + eps
  
  #put in data.frame for use in est.mu function
  dataY<-data.frame(time = t*dt, Y = y)
  
  mu[i] <- est.mu(dataY, hd, kern.leftexp, t.index = 500)$mu[1] #Take out a random timeindex (500)
}

mu_ <- sqrt(hd)*mu *sqrt(dt/2) #This should be N(0,ksq*omega^2) distributed from theorem 5

#Normalize
mu_normal <- sqrt(1/(ksq*sig2))*mu_

#test - it all looks really fine
mean(mu_normal)
var(mu_normal)
plot(mu_normal)
qqnorm(mu_normal)



############ Test convergence of sigma with noise (convergence of eq. 26 to equation in Theorem 6 p. 44) ############

#Simulate the BM with noise
x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
x <- cumsum(x)
eps <- rnorm(length(t), mean = 0, sd = sqrt(omega2))
y <- x + eps

#Put in dataframe
dataY<-data.frame(time = t*dt, Y = y)

#Get sigma estimates on 1000 different time_points (vigtigt at bandwidth er rigtig lille)
sigma_1000_times <- est.sigma(data = dataY, hv = hv, kern = kern.leftexp, wkern = kern.parzen,t.index = seq(1000,10000,by=10),lag=15)$sig

est_omega <- sqrt(dt)*sqrt(sigma_1000_times) #This should converge to K_2 * E(delta_eps^2) = K_2*2*E(eps^2) = K_2*2*omega2
est_one <- est_omega/(2*ksq*omega2) #This should converge in probability to 1

#Test that this is close to 1
plot(est_one)
mean(est_one)

#Right now the estimate just decreases when we increase the lag. 
# It looks as if there is missing some constant to counteract on this.

