setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")
source("Estimation/estimates.R")

tester<-function(testfun, norm = 1, sub = T, mode = FALSE, h = hd, t.index = 10000,
                 mean = 0, sd = sig, noise = 0,
                 t = 0:10000, mat = 1, N = 10000, plt = FALSE){
  # runs the testfun(ction) N times and returns mean/var/values/(plot)
  # sd is multiplied by sqrt(dt)
  
  dt <- mat/t[length(t)]
  out <- numeric(N)
  for(i in 1:N)
  {
    x <- rnorm(length(t) , mean = mean*dt, sd = sd*sqrt(dt))
    if(is.function(noise)){
      noise <- noise(t)
    }
    eps <- rnorm(length(t), mean = 0, sd = noise)
    x <- cumsum(x)
    y <- x + eps
    data<-data.frame(time = t*dt, Y = y)
    out[i]<-as.numeric(testfun(data, h, t.index = t.index, originalEstimator = mode)[2])
    out[i] <- norm*(out[i]-mean*sub)
    eps.last[i] <- eps[length(t)]
  }
  if(plt){
    qqnorm(out)
    abline(0,sd(out))
    
  }
  

  return(list(mean = mean(out), var = var(out), val = out, noise = noise, eps = eps.last))
}
# deets


mat <- 1
t <- 0:10000  # time_points
dt <- mat/t[length(t)]


sig <- 0.01
omega <- 0.5
ksq <- kern.leftexp$ksq # K2

hd <- 0.01 #bandwidth

C <- sqrt(ksq*omega^2) # K2/2 E (deps^2)

test4eps <- tester(est.mu, hd, T, T, mean = 0, noise = omega, plt = F)
test4eps$val
test4eps$eps

stop()
###########################
#       convergence       #
###########################

# NO drift | NO noise
#############################################
test1a <- tester(est.mu, sqrt(hd)/sqrt(sig^2*ksq), T, T, noise = 0, plt = T) # kim's
test1a$mean
test1a$var

#test1b <- tester(est.mu, sqrt(hd)/sqrt(sig^2*ksq), T, F, noise = 0, plt = T) # ours # WHAT TO NORMALIZE WITH?
#test1b$mean
#test1b$var
#############################################

# YES drift | NO noise
############################################# 
test2a <- tester(est.mu, sqrt(hd)/sqrt(sig^2*ksq), T, T, mean = 1, noise = 0, plt = T) # kim's
test2a$mean
test2a$var

#test2b <- tester(est.mu, sqrt(hd)/sqrt(sig^2*ksq), T, F, mean = 1, noise = 0, plt = T) # ours # WHAT TO NORMALIZE WITH?
#test2b$mean
#test2b$var
#############################################

# NO drift | YES noise
#############################################
test3a <- tester(est.mu, sqrt(hd*dt)/C, T, T, noise = omega, plt = T) # kim's 
test3a$mean #th5 p43
test3a$var

test3b <- tester(est.mu, sqrt(hd*dt)/C, T, F, noise = omega, plt = T) # ours
test3b$mean
test3b$var
#############################################

# YES drift | YES noise
#############################################
test4a <- tester(est.mu, sqrt(hd*dt)/C, T, T, mean = 1, noise = omega, plt = T) # kim's
test4a$mean
test4a$var

test4b <- tester(est.mu, sqrt(hd*dt)/C, T, F, mean = 1, noise = omega, plt = T) # ours
test4b$mean
test4b$var
#############################################

# EXTRA 
#############################################
test4a <- tester(est.mu, hd, T, T, mean = 0, noise = omega, plt = T) # kim's #different scaling
test4a$mean
test4a$var


test4b <- tester(est.mu, sqrt(hd*dt)/C, T, F, mean = 0, noise = omega*(dt*t), plt = T) # ours
test4b$mean
test4b$var
#############################################

# Sigma spot
#############################################
test5a <- tester(est.sigma, dt/(ksq*omega^2),F, T, noise = omega, plt = F) # kim's
test5a$mean
test5a$var

test5b <- tester(est.sigma, (dt/hd)/(C*sqrt(ksq*0.5)), F, F, noise = omega, plt = F) # uns
test5b$mean
test5b$var
#############################################