setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")
source("Estimation/estimates.R")

tester<-function(testfun, mode = FALSE, h = hd, t.index = 5000, mean = 0, sd = sig, noise = 0,
                 t = 0:10000, mat = 1, N = 1000, plt = FALSE){
  # runs the testfun(ction) N times and returns mean/var/values/(plot)
  # sd is multiplied by sqrt(dt)
  
  dt <- mat/t[length(t)]
  out = numeric(N)
  for(i in 1:N)
  {
    x <- rnorm(length(t) , mean = mean*dt, sd = sd*sqrt(dt))
    eps <- rnorm(length(t), mean = 0, sd = noise)
    x <- cumsum(x)
    y <- x + eps
    data<-data.frame(time = t*dt, Y = y)
    
    out[i]<-as.numeric(testfun(data, h, t.index = t.index, originalEstimator = mode)[2])
  
  }
  if(plt){
    qqnorm(out)
  }
  mv = c(mean(out), var(out))
  names(mv) = c("mean", "var")

  return(list(mv = mv))
}
# deets

mat <- 1
t <- 0:10000  # time_points
dt <- mat/t[length(t)]


sig <- 0.01
omega <- 0.5
ksq <- kern.leftexp$ksq # K2

hd <- 0.01 #bandwidth

Csq = ksq/2*omega # K2/2 E (deps^2)

# convergence
test1a <- tester(est.mu, T, noise = 0, plt = T) # kim's
sqrt(hd)*test1a$mv/(sig^2*ksq) # eqn 13 on p8

test1b <- tester(est.mu, F, noise = 0, plt = T) # ours
sqrt(hd*dt)*test1b$mv/Csq #noise is zero...

test2a <- tester(est.mu, T, noise = omega, plt = T) # kim's
sqrt(hd*dt)*test2a$mv/(ksq*omega^2) #th5 p43

test2b <- tester(est.mu, F, noise = omega, plt = T) # ours
sqrt(hd*dt)*test2b$mv/Csq

test3a <- tester(est.sigma, T, noise = omega, plt = T) # kim's
dt * test3a$mv/(ksq*omega^2) #p44

test3b <- tester(est.sigma, F, noise = omega, plt = T) # uns
dt/hd * test3b$mv/Csq

