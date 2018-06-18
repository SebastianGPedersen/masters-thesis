# TESTING AND BENCHMARKING
setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/pre-average.R")
source("kernels/kernels.R")
source("simulation/jumps.R")

# SIMULATION 
setting <- sim.setup(Npath = 2, Nstep = 23400, omega = 0.0000225)
hest <- sim.heston(setting)

sim <- sim.path(path = 1, sim.data = hest)
sim <- list(time = sim$time, Y = diff(sim$Y), raw = diff(sim$Y))

# ESTIMATION SETTINGS 
n = 23400
mat <- 6.5/(24*7*52)#*52*7*24*60*60 #In years
hd <- 300/(52*7*24*60*60)

# BENCHMARKERIA
require(ggplot2)
require(microbenchmark)
require(grid)
require(gridExtra)

N<-23

time <- matrix(NA, nrow = N, ncol = 5)

index <- 1000*(1:N)
for(i in 1:N){
  tind <- 1000*i
  if(i == 1){
    prev.mu <- est.mu.next(sim, hd = hd, t.index = tind)
    prev.sig <- est.sigma.next(sim, hv = hd, t.index = tind, lag = 10)
    time[i,] <- c(tind, summary(microbenchmark(est.mu(sim, hd = hd, t.index = tind),
                                       est.mu.next(sim, hd = hd, t.index = tind),
                                       est.sigma(sim, hv = hd, t.index = tind, lag = 10),
                                       est.sigma.next(sim, hv = hd, t.index = tind, lag = 10),
                                       times = 1000, unit = "ms"))$mean)
  }
  else{
    time[i,] <- c(tind, summary(microbenchmark(est.mu(sim, hd = hd, t.index = tind),
                                       est.mu.next(sim, prevmu = prev.mu, hd = hd, t.index = tind),
                                       est.sigma(sim, hv = hd, t.index = tind, lag = 10),
                                       est.sigma.next(sim, prevsig = prev.sig, hv = hd, t.index = tind, lag = 10),
                                       times = 1000, unit = "ms"))$mean)
    prev.mu <- est.mu.next(sim, hd = hd, t.index = tind)
    prev.sig <- est.sigma.next(sim, hv = hd, t.index = tind, lag = 10)
  }
}
data<-list(index = time[, 1])

data$index <- time[,1]; data$mu <- time[,2]; data$mu.optim <- time[,3]; data$sig <- time[,4]; data$sig.optim <- time[,5]

data<-data.frame(data)

g1 <- ggplot() +
  geom_line(data=data, aes(x=index, y=mu, col = "est.mu"), size = 1) +
  geom_line(data=data, aes(x=index, y=mu.optim, col = "est.mu.next"), size = 1) +
  ylab("Run time (ms)") +
  scale_x_continuous(name = "t.index", breaks = pretty(index, 2))

g2 <- ggplot() +
  geom_line(data=data, aes(x=index, y=sig, col = "est.sigma"), size = 1) +
  geom_line(data=data, aes(x=index, y=sig.optim, col = "est.sigma.next"), size = 1) +
  ylab("Run time (ms)") +
  scale_x_continuous(name = "t.index", breaks = pretty(index, 2))


grid.arrange(g1,g2,nrow = 1,
             top = textGrob("",gp=gpar(fontsize=20)))



# BENCHMARK THE OTHER METHODS #
# RCPP KNÆGTENE
#source("estimation/estimates_reloaded.R")
#source("estimation/estimates_revolution.R")

est.mu.mat.2.0 <- function(data, hd, kern = kern.leftexp, wkern=kern.parzen){
  
  #hd <- 5 / (60*24*7*52)
  
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hd)
  rescaling <- kern((data$time[2:length(data$time)]-t_now)/hd)
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  mus <- matrix(nrow = paths, ncol = n)
  
  dy <- data$Y
  
  for (path in 1:paths) {
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy
    
    #zero lag
    sum_terms <- products
    mu_non_scaled <- 1/hd * cumsum(sum_terms)
    mus[path,] <- mu_non_scaled/rescaling
  }
  
  return(list(time = data$time[-1],mu = mus)) #Don't include time zero, because dy doesn't include
}
est.sigma.mat.2.0 <- function(data, hv, kern = kern.leftexp, wkern=kern.parzen, lag=10){
  #data <- path
  #hv <- 5 / (60*24*7*52)
  #lag <- 10
  
  #p0 <- Sys.time()
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hv)
  rescaling <- kern((data$time[2:(length(data$time))]-t_now)/hv)
  
  #For sigma3.0
  #kernels <- kern((data$time[1:(length(data$time)-1)]-data$time[2:length(data$time)])/hv)
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  lags <- lag
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  #temp_func
  sigma_func <- function(path) {
    
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy
    
    # if (path == 1) {
    # print(dy[1:(lags+1)])
    # print(kernels[1:(lags+1)])
    # print(dy[1:(lags+1)]*kernels[1:(lags+1)])
    # print(sum((dy[1:(lags+1)]*kernels[1:(lags+1)])^2))
    # }
    
    #zero lag
    sum_terms <- products^2
    cum_sums <- cumsum(sum_terms)
    gamma_ls <- cum_sums
    sigmas <- gamma_ls[(lags+1):length(gamma_ls)] #It has to fit length-wise with all lags
    
    for (lag in 1:lags) {
      #lag <- 1
      sum_terms <- products[1:(n-lag)]*products[(lag+1):n]
      cum_sums <- cumsum(sum_terms)
      gamma_ls <- cum_sums
      
      sigmas <- sigmas + 2*wkern(lag,lags)*gamma_ls[(lags+1-lag):length(gamma_ls)]
    }
    sigmas <- 1/hv * sigmas/rescaling[(1+lags):n]^2 #rescaling
    
    return(sigmas)
  }
  
  sigmas <- foreach(path=1:paths, .combine = 'rbind')  %do% {sigma_func(path)}
  
  #Pad with zero in beginning for correct size
  if (paths == 1) {
    sigmas <- c(rep(0,lags),sigmas)
  } else {
    sigmas <- cbind(matrix(0,nrow = paths,ncol = lags),sigmas) 
  }
  
  #print(Sys.time()-p0)
  
  return(list(time = data$time[-1], sig = sigmas)) #Don't include time zero as this is not in dy
}

tind <- seq(10, 23400, by = 5)

# WAIT FOR 2.0 TO BE FIXED

# MU
microbenchmark(est.mu(sim, hd = hd, t.index = tind),
               est.mu.next(sim, hd = hd, t.index = tind),
               est.mu.mat.2.0(hest, hd),
               times = 100)

# SIGMA
microbenchmark(est.sigma(sim, hv = hd, t.index = tind, lag = 10),
               est.sigma.next(sim, hv = hd, t.index = tind, lag = 10),
               est.sigma.mat.2.0(hest, hd,lag = 10),
               est.sigma.next.cpp(sim, hv = hd, t.index = tind, lag = 10),
               times = 10)
