setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")
source("Estimation/estimates.R")

tester<-function(testfun, norm = 1, sub = T, mode = FALSE, h = hd, t.index = 10000,
                 mean = 0, sd = sig, noise = 0, rho=0,
                 t = 0:100000, mat = 1, N = 1000){
  # runs the testfun(ction) N times and returns mean/var/values/(plot)
  # sd is multiplied by sqrt(dt)
  
  simAR<-function(n, rho, sd){
    X<-rep(NA,n)
    X[1]<-rnorm(1, 0, sd/(1-rho))
    eps<-rnorm(n-1,0,sd)
    for(i in 2:n){
      X[i] <- rho*X[i-1]+eps[i-1]
    }
    return(X)
  }
  
  dt <- mat/t[length(t)]
  out <- numeric(N)
  eps.last <- numeric(N)
  for(i in 1:N)
  {
    x <- rnorm(length(t) , mean = mean*dt, sd = sd*sqrt(dt))
    x <- cumsum(x)
    if(rho != 0){
      eps <- simAR(length(t), rho, noise)
    }
    else{
      eps <- rnorm(length(t), mean = 0, sd = noise)
    }
    y <- x + eps
    data<-data.frame(time = t*dt, Y = y)
    out[i]<-as.numeric(testfun(data, h, t.index = t.index, originalEstimator = mode)[2])
    out[i] <- norm*(out[i]-mean*sub)
    eps.last[i] <- eps[t.index]
    if(i == 1){
      plot(y, type="l")
    }
  }
  if(plt){
    qqnorm(out)
    abline(0,sd(out))
    
  }
  return(list(mean = mean(out), var = var(out), val = out, noise = noise, eps = eps.last))
}

hd <- 0.01 #bandwidth

# JUST NOISE #
test0 <- tester(est.mu, hd, F, T, mean = 0, sd = 0.001, noise = 0.05, rho = 0, plt = F, last = 0.1)
test0$eps
test0$val
plot(test0$eps - test0$val)
plot(abs((test0$eps - test0$val)/test0$eps), ylim = c(0, 100))


# MILD NOISE #
test1 <- tester(est.mu, hd, F, T, mean = 1, sd = 1, noise = 0.01, rho = 0, plt = F, last = 0.1)
test1$eps
test1$val
plot(test1$eps - test1$val)
plot(abs((test1$eps - test1$val)/test1$eps), ylim = c(0, 100))


# MILD NOISE STRONG DIFF#
test2 <- tester(est.mu, hd, F, T, mean = 1, sd = 10, noise = 0.1, rho = 0, plt = F, last = 0.1)
test2$eps
test2$val
plot(test2$eps - test2$val)
plot((test2$eps - test2$val)/test2$eps, ylim = c(-100, 100))

# STRONG NOISE #
test3 <- tester(est.mu, hd, F, T, mean = 1, sd = 1, noise = 0.5, rho = 0, plt = F, last = 0.1)
test3$eps
test3$val
plot(test3$eps - test3$val)
plot((test3$eps - test3$val)/test3$eps, ylim = c(-100, 100))

# MILD AR #
test4 <- tester(est.mu, hd, F, T, mean = 1, sd = 1, noise = 0.01, rho = 0.1, plt = F, last = 0.1)
test4$eps
test4$val
plot(test4$eps - test4$val)
plot((test4$eps - test4$val)/test4$eps, ylim = c(-100, 100))

# STRONG AR #
test5 <- tester(est.mu, hd, F, T, mean = 1, sd = 1, noise = 0.01, rho = 0.5, plt = F, last = 0.1)
test5$eps
test5$val
plot(test5$eps - test5$val)
plot((test5$eps - test5$val)/test5$eps, ylim = c(-100, 100))
