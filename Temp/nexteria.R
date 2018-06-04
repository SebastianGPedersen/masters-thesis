setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("kernels/kernels.R")

require(Rcpp)

# START WITH MU
sourceCpp("estimation/next.cpp")


# faster versions
est.mu.next <-function(data, prevmu, hd, t.index){

    dy <- c(0,data$Y)
  
  kern <- kern.leftexp$kern
  # latest mu is picked out
  if(missing(prevmu)){
    prevmu <- list(time = -1, mu = 0)
  } 
  startmu<-list(time = prevmu$time[length(prevmu$time)], mu = prevmu$mu[length(prevmu$mu)])
  
  if(data$time[t.index][1] < startmu$time) stop("t.index should be higher than previous mu times")
  if(data$time[t.index][1] == startmu$time) t.index <- t.index[2:length(t.index)]
  
  t<-data$time[t.index]
  end <- t.index
  if(length(end) > 1){
    start <- c( min(which(data$time > startmu$time)), end[1:(length(end)-1)]+1)
  }
  else{
    start <- min(which(data$time > startmu$time))
  }
  
  # scaling
  dt1<-t[1]-startmu$time
  dt<-c(dt1, diff(t))
  rescale <- exp(-dt/hd)
  
  n <- length(data$time)
  
  # init
  tt = length(t)
  mu <- numeric(tt)
  # CALCULATE NEXT 'BLOCK'
  
  
  mu[1] <- startmu$mu*rescale[1] + 1/hd * sum( kern( (data$time[start[1]:end[1]] - t[1])/hd)*dy[(start[1]+1):(end[1]+1)] )
  if(tt > 1){
    for(j in 2:(tt)){
      mu[j] <- mu[j-1]*rescale[j] + 1/hd * sum( kern(  (data$time[start[j]:end[j]] - t[j])/hd)*dy[(start[j]+1):(end[j]+1)] )
    }
  }
  return(list(time = t, mu = mu))
}

# faster versions
est.mu.next.cpp <-function(data, prevmu, hd, t.index){
  
  dy <- c(0,data$Y)
  
  kern <- kern.leftexp$kern
  # latest mu is picked out
  if(missing(prevmu)){
    prevmu <- list(time = -1, mu = 0)
  } 
  startmu<-list(time = prevmu$time[length(prevmu$time)], mu = prevmu$mu[length(prevmu$mu)])
  
  if(data$time[t.index][1] < startmu$time) stop("t.index should be higher than previous mu times")
  if(data$time[t.index][1] == startmu$time) t.index <- t.index[2:length(t.index)]
  
  t<-data$time[t.index]
  end <- t.index
  if(length(end) > 1){
    start <- c( min(which(data$time > startmu$time)), end[1:(length(end)-1)]+1)
  }
  else{
    start <- min(which(data$time > startmu$time))
  }
  
  # scaling
  dt1<-t[1]-startmu$time
  dt<-c(dt1, diff(t))
  rescale <- exp(-dt/hd)
  
  
  mu <- mu_internal_cpp(time = data$time, dy = dy, t = t, rescale = rescale,
                        start = start, end = end, hd = hd, startmu = startmu$mu)
  
  
  return(list(time = t, mu = mu))
}

est.sigma.next.cpp <- function(data, prevsig, hv, t.index, wkern=kern.parzen, lag="auto"){   #
  # data list should include a times column and the Y column (log returns)
  # Handle lag
  if(lag=="auto") lag = 15 #temp
  
  if(t.index[1] < lag) stop("t.index can and should NOT be lower than lag length!")
  
  kern <- kern.leftexp$kern
  wkern = kern.parzen$kern
  
  dy <- c(0,data$Y)
  
  miss = FALSE
  if(missing(prevsig)){
    # initial handling is pretty weird because of lag length
    prevsig <- est.sigma(data, hv=hv, t.index = t.index[1],lag = lag) # change to minus
    miss = TRUE
    if(length(t.index) == 1){
      return(prevsig)
    }
  } 
  startsig<-list(time = prevsig$time[length(prevsig$time)], sig = prevsig$sig[length(prevsig$sig)])
  
  if(data$time[t.index][1] < startsig$time) stop("t.index should be higher than previous mu times")
  if(data$time[t.index][1] == startsig$time) t.index <- t.index[2:length(t.index)]
  
  t<-data$time[t.index]
  end <- t.index
  if(length(end) > 1){
    start <- c( min(which(data$time > startsig$time)), end[1:(length(end)-1)]+1)
  }
  else{
    start <- min(which(data$time > startsig$time))
  }
  
  # scaling
  dt1<-t[1]-startsig$time
  dt<-c(dt1, diff(t))
  rescale <- exp(-2*dt/hv)
  
  sig = sig_internal_cpp(data$time, dy, t, rescale, start, end, hv, lag, startsig$sig)
  
  if(miss){
    t <- c(prevsig$time, t)
    sig <- c(prevsig$sig, sig)
  }
  
  # return list
  return(list(time = t, sig = sig))
}

data <- list(time = 1:100*0.01, Y = 1:99)
hd <- 0.7
t.index <- c(5,10,99)

est.mu(data = data, hd = hd, t.index = t.index)

est.mu.next(data, hd = hd, t.index = t.index)

est.mu.next.cpp(data, hd = hd, t.index = t.index)

# SIGMARIA

prev = list(time = 0, sig = 0)

est.sigma(data, hd, t.index, lag = 3)


est.sigma.next.cpp(data, prev, hv = hd, t.index = t.index, lag = 3)

