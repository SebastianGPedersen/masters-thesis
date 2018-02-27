setwd(Sys.getenv("masters-thesis"))

source("estimation/estimates.R")

# Time setup
t <- 0:1000  # time
dt <- 1/t[length(t)]

# Define sigma
sig2 <- 1
omega = 10
ksq = 0.5 # K2

# omega = ((1/2)*sqrt(sig2*dt/t[length(t)]))^2


# first, simulate a set of random deviates
y <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt)) + rnorm(length(t), mean = 0, sd = omega)
y <- cumsum(y)

# form data and plot
plot(y)
dataX<-data.frame(time = t, Y = y)

# estimation of sigma
sig<-est.sigma(data = dataX, hv = 0.02, kern = kern.leftexp, wkern = kern.parzen, t.index = 500)

# dt*sigma should converge in P to omega^2
dt*sig$sig
ksq*omega^2


# test of convergence of mu to N(0, ksq*omega^2)
{
  hd<-0.02
  N = 1000
  mu = numeric(N)
  for(i in 1:N)
  {
    sig2 <- 1
    
    y <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt)) + rnorm(length(t), mean = 0, sd = omega)
    x <- cumsum(y)
    
    dataX<-data.frame(time = t*dt, Y = x)
    
    mu[i] <- sqrt(dt)*sqrt(hd)*(est.mu(dataX, hd, kern.leftexp, t.index = 500)$mu[1])
  }
  # properties of mu  - (   N(0,ksq*omega^2)   )
  mean(mu)
  var(mu)
  plot(mu)
  qqnorm(mu)
}

