setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("kernels/kernels.R")

# needs something to handle every second - right now it is assumed that data can be used raw
# last point is abit weird...
est.mu.next <-function(data, prevmu, hd, t.index){
  # data should contain time | logreg - prevmu should contain time | mu
  # t.index indicates where we wish to calculate the estimator
  
  # --- debug ---
  # data<-list(time = 1:10/10, Y = 1:10)
  # hd = 0.1
  # prevmu <- est.mu(data, hd, t.index = 1:3, originalEstimator = T)
  # t.index<-c(3,5,9)
  
  
  kern <- kern.leftexp$kern
  # latest mu is picked out
  if(missing(prevmu)){
    # neg time such that it is guarentee that it comes before any data that could start at t = 0
    # mu is zero so any weird time scaling will be irrelevant anyway
    prevmu <- list(time = -1, mu = 0)
  } 
  else{
    startmu<-list(time = prevmu$time[length(prevmu$time)], mu = prevmu$mu[length(prevmu$mu)])
  }
  
  # CALCULATE NEXT 'BLOCK'
  #nextdata <- list(time = data$time[data$time>startmu$time], Y = data$Y[data$time>startmu$time]) #>= or > hm?
  
  if(data$time[t.index][1] < startmu$time) stop("t.index should be higher than previous mu times")
  if(data$time[t.index][1] == startmu$time) t.index <- t.index[2:length(t.index)]
  
  t<-data$time[t.index]
  end <- t.index
  start <- c( min(which(data$time > startmu$time)), end[1:(length(end)-1)]+1)
  
  # scaling
  dt1<-t[1]-startmu$time
  dt<-c(dt1, diff(t))
  rescale <- exp(-dt/hd)
  
  
  dy <- diff(data$Y)
  if(anyNA(dy[t.index])) warning("t.index cannot be higher than t_n-1 (this will rightfully give you NAs)")
  
  n <- length(data$time)
  
  # init
  tt = length(t)
  mu <- numeric(tt)
  
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
a

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

setting <- sim.setup(Npath = 2, Nstep = 10000)
sims<-sim.heston(setting)
sim<-sim.path(1, sims)

(tind<-seq(from=1000, to=8000, by=100))

(test<-est.mu(data = sim, hd = 0.001, kern = kern.leftexp, t.index = tind,  originalEstimator = T)$mu
  -est.mu.next(data = sim, hd = 0.001, t.index = tind)$mu)

require(microbenchmark)

microbenchmark(est.mu(data = sim, hd = 0.001, t.index = tind, originalEstimator = T),
               est.mu2(data = sim, hd = 0.001, t.index = tind, originalEstimator = T),
               est.mu.next(data = sim, hd = 0.001, t.index = tind), times = 100)
