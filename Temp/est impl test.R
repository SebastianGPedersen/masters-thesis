setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("kernels/kernels.R")

est.mu.next <-function(data, prevmu, hd, t.index, originalEstimator = F){
  # data should contain time | logreg - prevmu should contain time | mu
  # t.index indicates where we wish to calculate the estimator
  
  # --- debug ---
  # data<-list(time = 1:10/10, Y = 1:10)
  # hd = 0.1
  # prevmu <- est.mu(data, hd, t.index = 1:3, originalEstimator = T)
  # t.index<-c(3,5,9)
  
  if(is.null(data$dy)){
    dy <- diff(data$Y)
  } else {
    dy <- data$dy
  }
  
  kern <- kern.leftexp$kern
  # latest mu is picked out
  if(missing(prevmu)){
    # neg time such that it is guarentee that it comes before any data that could start at t = 0
    # mu is zero so any weird time scaling will be irrelevant anyway
    prevmu <- list(time = -1, mu = 0)
  } 
  startmu<-list(time = prevmu$time[length(prevmu$time)], mu = prevmu$mu[length(prevmu$mu)])
  
  if(data$time[t.index][1] < startmu$time) stop("t.index should be higher than previous mu times")
  if(data$time[t.index][1] == startmu$time) t.index <- t.index[2:length(t.index)]
  
  t<-data$time[t.index]
  end <- t.index
  start <- c( min(which(data$time > startmu$time)), end[1:(length(end)-1)]+1)
  
  # scaling
  dt1<-t[1]-startmu$time
  dt<-c(dt1, diff(t))
  rescale <- exp(-dt/hd)
  
  if(anyNA(dy[t.index])) warning("t.index cannot be higher than t_n-1 (this will rightfully give you NAs)")
  
  n <- length(data$time)
  
  # init
  tt = length(t)
  mu <- numeric(tt)
  # CALCULATE NEXT 'BLOCK'
  mu[1] <- startmu$mu*rescale[1] + 1/hd * sum( kern( (data$time[start[1]:end[1]] - t[1])/hd)*dy[start[1]:end[1]] )
  for(j in 2:(tt)){
    mu[j] <- mu[j-1]*rescale[j] + 1/hd * sum( kern(  (data$time[start[j]:end[j]] - t[j])/hd)*dy[start[j]:end[j]] )
  }
  
  # debug checker
  #1/hd * sum( kern( (data$time[1:(n-1)] - data$time[5])/hd)*dy[1:(n-1)] )
  
  #ADD BLOCK TO EXISTING (wow much blockchainy) (if it existed)
  if(min(prevmu$time) >= 0){
    t <- c(prevmu$time, t)
    mu <- c(prevmu$mu, mu)
  } 
  
  # return list
  return(list(time = t, mu = mu))
}

est.sigma.next <- function(data, prevsig, hv, t.index, wkern=kern.parzen, lag="auto"){   #
  # data list should include a times column and the Y column (log returns)
  
  # Handle lag
  if(lag=="auto") lag = 15 #temp
  
  if(t.index[1] < lag) stop("t.index can and should NOT be lower than lag length!")
  
  kern <- kern.leftexp$kern
  wkern = kern.parzen$kern
  
  if(is.null(data$dy)){
    dy <- diff(data$Y)
  } else {
    dy <- data$dy
  }
  
  # --- debug ---
  {
    #data<-list(time = 1:100/100, Y = 1:100)
    #hv = 0.1
    #prevsig <- est.sigma(data, hv, t.index = 2:4, originalEstimator = T, lag = lag)
    #t.index<-c(5,9)
  }
  
  if(missing(prevsig)){
    # initial handling is pretty weird because of lag length
    prevsig <- est.sigma(data, hv, t.index = t.index[1], originalEstimator = T, lag = lag) # change to minus
  } 
  startsig<-list(time = prevsig$time[length(prevsig$time)], sig = prevsig$sig[length(prevsig$sig)])
  
  if(data$time[t.index][1] < startsig$time) stop("t.index should be higher than previous mu times")
  if(data$time[t.index][1] == startsig$time) t.index <- t.index[2:length(t.index)]
  
  t<-data$time[t.index]
  end <- t.index
  start <- c( min(which(data$time > startsig$time)), end[1:(length(end)-1)]+1)
  
  # scaling
  dt1<-t[1]-startsig$time
  dt<-c(dt1, diff(t))
  rescale <- exp(-2*dt/hv)
  
  if(anyNA(dy[t.index])) warning("t.index cannot be higher than t_n-1 (this will rightfully give you NAs)")
  
  n <- length(data$time)
  
  # init
  tt = length(t)
  sig <- numeric(tt)
  
  #Define gamma
  gamma<-function(l){
    out<-sum(   kern( (data$time[start[j]:end[j]] - t[j])/hv )    *dy[start[j]:end[j]]*       
                  kern( (data$time[(start[j]-l):(end[j]-l)] - t[j])/hv )*dy[(start[j]-l):(end[j]-l)]   )
    return(out)
  }
  # CALCULATE NEXT 'BLOCK'
  j<-1
  sig[1] <- startsig$sig*rescale[1] + sum(  (kern(   (data$time[start[1]:end[1]] - t[1])/hv   )*dy[start[1]:end[1]])^2 )  # l = 0
  if (lag >=1) {
    for(l in 1:lag){
      sig[1] <- sig[1] + 2*(wkern(l/(lag+1))*gamma(l) )
    }
  }
  for(j in 2:tt){
    sig[j] <- sig[j-1]*rescale[j] + sum(  (kern(   (data$time[start[j]:end[j]] - t[j])/hv   )*dy[start[j]:end[j]])^2  )  # l = 0
    if (lag >=1) {
      for(l in 1:lag){
        sig[j] <- sig[j] + 2*(wkern(l/(lag+1))*gamma(l) )
      }
    }
  }
  
  # --- debug ---
  {
    # husk at j bruges i udregning 1 - husk at k�r j<-1 hver gang der tjekkes op mod noget i sig[1]!!!!
    
    #sum( (kern( (data$time[1:(n-1)] - t[1])/hv)*dy[1:(n-1)])^2 )+2*wkern(1/2)*gamma2(1, t[1])
    #startsig$sig<-sum( (kern( (data$time[1:(n-1)] - data$time[4])/hv)*dy[1:(n-1)])^2 )+2*wkern(1/2)*gamma2(1, data$time[4])
    
    #gamma2<-function(l, t){
    #  out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )*dy[(l+1):(n-1)]*       
    #                 kern( (data$time[1:(n-l-1)] - t)/hv )*dy[1:(n-l-1)]   )
    #  return(out)
    #}
  }
  
  #ADD BLOCK TO EXISTING (wow much blockchainy) (if it existed)
  t <- c(prevsig$time, t)
  sig <- c(prevsig$sig, sig)
  
  # return list
  return(list(time = t, sig = sig))
}


est.mu2 <- function(data, hd, kern = kern.leftexp, t.index, t.points, originalEstimator=FALSE){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  n = length(data$time)
  
  # mode-handling
  mode = 0
  if(missing(t.index) & missing(t.points)){
    mode = 1
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind<-1:length(data$time)
  }
  else if(missing(t.index) & !missing(t.points)){
    mode = 2
    t<-t.points
    ind = numeric(length(t))
    #ind[1] = sum(data$time < t[1])
    #ind[i] = data$time[ind[i-1]+sum(data$time[ind[i-1]:n] < t[i] )] #for 2 to n this may be faster than which.max
    for(i in 2:length(t)){
      ind[i] = which.max(data$time[data$time<=t[i]])
    }
  }
  else{
    mode = 3
    t<-data$time[t.index]
    #if(min(t.index) < 2){
    #stop("t.index cannot be less than two (deltaY[1] is not defined)")
    #}
    ind<-t.index
  }
  # t should now be data$time points
  
  if(!originalEstimator){
    diffY<-diff(data$Y)
    everySecondSeq<-seq(1L, length(diffY), by=2L)
    dy<-diffY[everySecondSeq]
    tempTime <- data$time[seq(1L, length(data$time), by = 2L)] #temporary placeholder, for compatibility
    data <- NULL #remove data. Handles errors with e.g. data.frames
    data$time <- tempTime #creates list compatable with below
    
    # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  } else {
    
    dy = diff(data$Y)
  }
  #update n accordingly
  n = length(data$time)
  tt = length(t)
  mu = numeric(tt)
  
  for(j in 1:(tt-1)){
    mu[j] = (1/hd)*sum(kern((data$time[1:(ind[j])] - t[j])/hd)*dy[1:ind[j]])
  }
  mu[tt] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[tt])/hd)*dy[1:(n-1)])
  
  # return list
  return(list(time = t, mu = mu))
}
est.sigma2 <- function(data, hv, kern, wkern, t.index, lag="auto"){   # we could do lag = "auto"
  # data list should include a times column and the Y column (log returns)
  
  # Handle lag
  if(lag=="auto") lag = 15 #temp
  
  t<-data$time[t.index]
  ind<-t.index
  
  tt = length(t)
  sig = numeric(tt) 
  
  
  dy = diff(data$Y)
  
  gamma2<-function(l, t){
    # to be used in the loop with j as n            
    out<-ifelse(ind[j]+l < n,
                sum(   kern( (data$time[(l+1):(ind[j]+l)] - t)/hv )*         #if true
                         dy[(l+1):(ind[j]+l)]*       
                         kern( (data$time[1:ind[j]] - t)/hv )*
                         dy[1:ind[j]]   ),
                sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )*              #if false
                         dy[(l+1):(n-1)]*       
                         kern( (data$time[1:(n-l-1)] - t)/hv )*
                         dy[1:(n-l-1)]   ))
    return(out)
  }
  
  gamma<-function(l, t){
    # the standard where end = n
    l <- abs(l)
    out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )* 
                   dy[(l+1):(n-1)]*       #indices are literally #1 reason for bugs
                   kern( (data$time[1:(n-l-1)] - t)/hv )*
                   dy[1:(n-l-1)]   )
    return(out)
  }
  
  for (j in 1:(tt-1)) {
    sig[j] = sum(  (kern(   (data$time[1:ind[j]] - t[j])/hv   )*dy[1:ind[j]])^2  )  # l = 0
    if (lag >=1) {
      for(l in 1:lag){
        sig[j] = sig[j] + 2*(wkern(l/(lag+1))*gamma2(l,t[j])) #gamma2
        
      }
    }
  }
  sig[tt] = sum(  (kern(   (data$time[1:(n-1)] - t[tt])/hv   )*dy[1:(n-1)])^2  )  # l = 0
  if (lag >=1) {
    for(l in 1:lag){
      sig[tt] = sig[tt] + 2*(wkern(l/(lag+1))*gamma(l, t[tt]))
    }
  }
  
  # return list
  return(list(time = t, sig = sig))
}

setting <- sim.setup(Npath = 2, Nstep = 23400)
sims<-sim.heston(setting)
sim<-sim.path(1, sims)

plot(sims$vol[1,], type = "l")

(tind<-seq(from=60, to=10000, by=60))
# ~~~~~~~~~~~~~~~~~~ MU ~~~~~~~~~~~~~~~~

data2 <- est.EveryOtherDiffData(sim)

(test<-est.mu(data = data2, hd = 0.0001, kern = kern.leftexp, t.index = tind,  originalEstimator = T)$mu
  -est.mu.next(data = data2, hd = 0.0001, t.index = tind, originalEstimator = T)$mu)

(timer<-est.mu(data = data2, hd = 0.0001, kern = kern.leftexp, t.index = tind,  originalEstimator = T)$time
  -est.mu.next(data = data2, hd = 0.0001, t.index = tind, originalEstimator = T)$time)


# ~~~~~~~~~~~~~~~~~~SIGMA ~~~~~~~~~~~~~~
(lags<-(est.sigma(lag=0, data = data2, hv = 0.0001, t.index = tind, originalEstimator = T)$sig
        -est.sigma.next(lag=0, data = data2, hv = 0.0001, t.index = tind)$sig))

(timematch<-est.sigma(lag=0, data = sim, hv = 0.0001, t.index = tind, originalEstimator = T)$time
  -est.sigma.next(lag=0, data = sim, hv = 0.0001, t.index = tind)$time)

plot(est.sigma(lag=0, data = sim, hv = 0.0001, t.index = tind, originalEstimator = T)$sig, type = "l")
lines(est.sigma.next(lag=0, data = sim, hv = 0.0001, t.index = tind)$sig, col = "red")

(netto<-(est.sigma(data = data2, hv = 0.0001, kern = kern.leftexp, wkern = kern.parzen, t.index = tind, originalEstimator = T)$sig
         -est.sigma.next(data = data2, hv = 0.0001, t.index = tind)$sig))

require(microbenchmark)

(comp.mu<-microbenchmark(est.mu(data = sim, hd = 0.01, t.index = tind, originalEstimator = T),
                         est.mu2(data = sim, hd = 0.01, t.index = tind, originalEstimator = T),
                         est.mu.next(data = sim, hd = 0.01, t.index = tind), times = 100))

(comp.sig<-microbenchmark(est.sigma(data = sim, hv = 0.01, t.index = tind, originalEstimator = T)$sig,
                          -est.sigma.next(data = sim, hv = 0.01, t.index = tind)$sig, times = 5))
