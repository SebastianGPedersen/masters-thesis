# HESTON SIMULATION PLAYGROUND

setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")

sim.heston.reflect<-function(settings){
  
  N = settings$Npath
  mat = settings$mat
  steps = settings$Nsteps
  kappa = settings$kappa
  c = settings$c
  theta = settings$theta
  xi = settings$xi       ##vol of vol
  rho = settings$rho
  omega = settings$omega
  
  dt = mat/steps #dt is in years
  time = 0:steps*dt # maybe time should be relative (0-1)
  #time = 0:(steps)/(steps)
  
  X = matrix(nrow = N, ncol = steps+1)
  Y = matrix(nrow = N, ncol = steps+1)
  vol = matrix(nrow = N, ncol = steps+1) #matrix if we want to save values along the way
  call = numeric(N) # testing purpose
  
  X[, 1] = 0
  vol[, 1] = rgamma(N, 2*kappa*theta/xi^2, 2*kappa/xi^2)
  #Y[, 1] = X[,1] + gamma*sqrt(vol[,1])/sqrt(steps)*rnorm(N,0,1) #Changed from vol to sqrt(vol) /Seb 20.02.18
  Y[ ,1] = X[, 1] + omega*rnorm(N,0,1)
  
  for(i in 2:(steps+1)){
    NS = rnorm(N,0,1)
    NV = rho*NS + sqrt(1-rho^2)*rnorm(N,0,1) #From StatÃ (Olivier) Theorem I.5 or Graphical example 1.20 

    X[,i] =   X[,i-1] + c*vol[,i-1]*dt + sqrt(vol[,i-1])*sqrt(dt)*NS                          # logs (driftless)
    vol[, i] = abs(vol[,i-1]+kappa*(theta-vol[,i-1])*dt + xi*sqrt(vol[,i-1]*dt)*NV)
    
    #Observed Y
    #omega = gamma*sqrt(vol[,i])/sqrt(steps)     # n corresponds to steps and not repetitions N? #should vol be i-1? No?
    Y[,i] = X[,i] + omega * rnorm(N,0,1)
    
  }
  return(list(time = time, Y = Y, X = X, vol = vol))
}

setting <- sim.setup(Npath = 1)

set.seed(2342)
a<-sim.heston(setting)
set.seed(2342)
b<-sim.heston.reflect(setting)

plot(a$X[1,], type = "l")
lines(b$X[1,], type = "l", col = "red")

plot(a$vol[1,], type = "l")
lines(b$vol[1,], type = "l", col = "red")
