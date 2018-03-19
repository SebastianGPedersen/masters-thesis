setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")
#source("Estimation/estimates.R")

# Time parameters
mat <- 1
t <- 0:10000  # time_points
dt <- mat/t[length(t)]

# Define sigma
sig <- 0.01
sig2 <- sig^2
omega <- 0.5
omega2 <- omega^2
ksq <- 0.5 # K2
hd <- 0.01 #bandwidth in mu
hv <- 0.01 #bandwidth in sigma


# Simulate a Brownian motion (X) and one with microstrucutre noise (Y) ~ N(0, omega^2)
x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
eps <- rnorm(length(t), mean = 0, sd = sqrt(omega))
x <- cumsum(x)
y <- x + eps


# t<-1:10
# y<-2^(1:10)
# dataY<-data.frame(time=t,Y=y)
# 
# dataY
# est.mu(dataY, hd, kern.leftexp, t.index=8, t.points=NA, originalEstimator=F)
# Plot part of x and y to validate it is done correctly
plot(x[1:100],type="l")
plot(y[1:100],type = "l")
#Støjen overgår tydeligt X som den burde


###################################################################################################################
###################################### Test convergence of mu with noise  #########################################
###################################################################################################################

N <- 1000
mu <- numeric(N)
sigmaEst <- numeric(N)


######################################### Normal noise ######################################
simExpr<-substitute(rnorm(length(t), mean = 0, sd = sqrt(omega2)))

for(i in 1:N) #Takes approx 30sec
{
  #simulate brownian motion with noise
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x)
  eps <- eval(simExpr)
  y <- x + eps
  
  #put in data.frame for use in est.mu function
  dataY<-data.frame(time = t*dt, Y = y)
  
  mu[i] <- est.mu(dataY, hd, kern.leftexp, t.index = 9000)$mu[1] #Take out a random timeindex (500)
}

mu_ <- sqrt(dt*hd)*mu #This should be N(0, (K2/2)*(2*omega2) )

#Normalize
mu_normal <- sqrt((1/((K2/2)*(2*omega2))))*mu_

#test - it all looks really fine
mean(mu_normal)
var(mu_normal)
plot(mu_normal)
qqnorm(mu_normal)
qqline(mu_normal)

######################################### t-distributed noise ######################################

degreesOfFreedom <- 3
simExpr<-substitute(rt(length(t), df = degreesOfFreedom))

for(i in 1:N) #Takes approx 30sec
{
  #simulate brownian motion with noise
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x)
  eps <- eval(simExpr)
  y <- x + eps
  
  #put in data.frame for use in est.mu function
  dataY<-data.frame(time = t*dt, Y = y)
  
  mu[i] <- est.mu(dataY, hd, kern.leftexp, t.index = 9000)$mu[1] #Take out a random timeindex (9000)
}

mu_ <- sqrt(dt*hd)*mu #This should be N(0, (K2/2)*2*dof/(dof-1)) 
tVar<-degreesOfFreedom/(degreesOfFreedom-2)
#Normalize
mu_normal <- sqrt((1/((K2/2)*(2*tVar))))*mu_

#QQ-plot
mean(mu_normal)
var(mu_normal)
plot(mu_normal)
qqnorm(mu_normal)
qqline(mu_normal)


######################################################################################################################
###################################### Test convergence of sigma with noise  #########################################
######################################################################################################################

######################################### Normal noise ######################################
simExpr<-substitute(rnorm(length(t), mean = 0, sd = sqrt(omega2)))

for(i in 1:N) #Takes approx 30sec
{
  #simulate brownian motion with noise
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x)
  eps <- eval(simExpr)
  y <- x + eps
  
  #put in data.frame for use in est.mu function
  dataY<-data.frame(time = t*dt, Y = y)
  
  #Get sigma estimates (vigtigt at bandwidth er rigtig lille)
  sigmaEst[i] <- est.sigma(data = dataY, hv = hv, kern = kern.leftexp, wkern = kern.parzen,t.index = 9000,lag=15)$sig
  
}

est_omega <- (dt/hv)*sigmaEst #This should converge to (K_2/2) * E(delta_eps^2) = (K_2/2)*2*E(eps^2) = (K_2/2)*2*omega2
est_one <- est_omega/((K2/2)*2*omega2) #This should converge in probability to 1
# plot(est_omega)
# mean(est_omega)
# var(est_omega)

#Test that this is close to 1
plot(est_one)
mean(est_one)
var(est_one)

######################################### t-distributed noise ######################################
degreesOfFreedom <- 3
simExpr<-substitute(rt(length(t), df = degreesOfFreedom))

for(i in 1:N) #Takes approx 30sec
{
  #simulate brownian motion with noise
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x)
  eps <- eval(simExpr)
  y <- x + eps
  
  #put in data.frame for use in est.mu function
  dataY<-data.frame(time = t*dt, Y = y)
  
  #Get sigma estimates (vigtigt at bandwidth er rigtig lille)
  sigmaEst[i] <- est.sigma(data = dataY, hv = hv, kern = kern.leftexp, wkern = kern.parzen,t.index = 9000,lag=15)$sig
  
}
mu_ <- sqrt(dt*hd)*mu

est_omega <- (dt/hv)*sigmaEst #This should converge to (K_2/2) * E(delta_eps^2) = (K_2/2)*2*E(eps^2) = (K_2/2)*2*(K2/2)*2*dof/(dof-1)
#Normalize
tVar<-degreesOfFreedom/(degreesOfFreedom-2)
est_one <- est_omega/((K2/2)*2*tVar) #This should converge in probability to 1
# plot(est_omega)
# mean(est_omega)
# var(est_omega)

#Test that this is close to 1
plot(est_one)
mean(est_one)
var(est_one)

######################################### Normal Correlated noise ######################################
simAR<-function(n, rho, sd){
  X<-rep(NA,n)
  X[1]<-rnorm(1, 0, sd/(1-rho))
  eps<-rnorm(n-1,0,sd)
  for(i in 2:n){
    X[i] <- rho*X[i-1]+eps[i-1]
  }
  return(X)
}

rho<-0.1
for(i in 1:N) #Takes approx 30sec (i.e. slow)
{
  #simulate brownian motion with noise
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x)
  eps <- simAR(length(t), rho, sqrt(omega2))
  y <- x + eps
  
  #put in data.frame for use in est.mu function
  dataY<-data.frame(time = t*dt, Y = y)
  
  #Get sigma estimates (vigtigt at bandwidth er rigtig lille)
  mu[i] <- est.mu(dataY, hd, kern.leftexp, t.index = 9000)$mu[1] #Take out a random timeindex (9000)
  sigmaEst[i] <- est.sigma(data = dataY, hv = hv, kern = kern.leftexp, wkern = kern.parzen,t.index = 9000,lag=15)$sig
  
}

estT<-hd*(mu/sqrt(sigmaEst))


mean(estT)
var(estT)
plot(estT)
qqnorm(estT)
qqline(estT)


# estT2<-hd*(mu/sigmaEst)
# mean(estT2)
# var(estT2)
# plot(estT2)
# qqnorm(estT2)
# qqline(estT2)
# 









mu_ <- sqrt(dt*hd)*mu
est_omega <- (dt/hv)*sigmaEst #Thid converge to (K_2/2) * E(delta_eps^2) 
est_one <- est_omega/((K2/2)*(1+(1-rho)^2+(1-rho)^2*rho^2/(1-rho^2))*omega2) #This should converge in probability to 1
#((K2/2)*(1+(1-rho)^2+(1-rho)^2*rho^2/(1-rho^2))*omega2)
# plot(est_omega)
#mean(est_omega)
# var(est_omega)

#Test that this is close to 1
plot(est_one)
mean(est_one)
var(est_one)

