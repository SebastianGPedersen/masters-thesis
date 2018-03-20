
# BANDWIDTH SHOULD BE TRANSLATED FROM SECONDS TO YEARS ( BW / Seconds per year)
# Consider multiplying time such that our unit is in seconds and not years...!

est.mu <- function(data, hd, kern = kern.leftexp, t.index, t.points, originalEstimator=FALSE){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # t.points does NOT work yet!
  
  if(is.null(data$dy)){
    dy <- diff(data$Y)
  } else {
    dy <- data$dy
  }
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  # mode-handling
  mode = NA
  if(missing(t.index) & missing(t.points)){
    mode = 1
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind = 1:(length(data$time))
  }
  else if(missing(t.index) & !missing(t.points)){
    mode = 2
    t<-t.points
    ind = numeric(length(t))
    #ind[1] = sum(data$time < t[1])
    #ind[i] = data$time[ind[i-1]+sum(data$time[ind[i-1]:n] < t[i] )] #for 2 to n this may be faster than which.max
    for(i in 2:length(t)){
      ind[i] = which.max(data$time[data$time<t[i]])
    }
  }
  else{
    mode = 3
    t<-data$time[t.index]
    ind = t.index
  }
  
  #t should now be data$time points
  #n = length(t)
  tt = length(t)
  n = length(data$time)
  mu = numeric(tt)          # We can only have bandwidth to end amount of calcs
  
  # Optimization removed
  if(mode == 1){
    for(j in 1:tt){
      mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[1:(n-1)])   
    }
  }
  else if(mode == 3){
    for(j in 1:tt){
      mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[1:(n-1)])
    }
  }
  else{ # time points not implemented
    mu[j] = NA
  }
  
  # return list
  return(list(time = t, mu = mu))
}

est.sigma <- function(data, hv, kern = kern.leftexp, wkern = kern.parzen, t.index, lag="auto", originalEstimator=FALSE){   # we could do lag = "auto"
  # data list should include a times column and the Y column (log returns)
  
  # Handle lag
  if(lag=="auto") lag = 15 #temp
  
  if(is.null(data$dy)){
    dy <- diff(data$Y)
  } else {
    dy <- data$dy
  }
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  # if missing handling:
  if(missing(t.index)){
    start = lag+1
    t<-data$time[start:(length(data$time))] # if nothing specified - every point in data
    ind<-start:(length(data$time))
  }else{
    t<-data$time[t.index]
    ind<-t.index
  }
  # t should now be data$time points

  n = length(data$time)
  tt = length(t)
  sig = numeric(tt) 
  
  gamma<-function(l, t, end){
    # end is the highest needed index
    out<- sum(   kern( (data$time[(l+1):(end-1)] - t)/hv )* 
                  dy[(1+l):(end-1)]*       
                 kern( (data$time[1:(end-l-1)] - t)/hv )*
                    dy[1:(end-1-l)]   )
    return(out)
  }
  

  end = n
  for (j in 1:tt) {
   sig[j] = sum(  (kern(   (data$time[1:(end-1)] - t[j])/hv   )*
                   dy[1:(end-1)])^2  )  # l = 0
   
   #sig[j] = sum (  dy[2:end]^2  )
   if (lag >=1) {
     for(l in 1:lag){
       #sig[j] = sig[j] + 2*(wkern(l/n)*gamma(l,t[j], end))
       sig[j] = sig[j] + 2*(wkern(l/(lag+1))*gamma(l,t[j], end))
       #sig[j] = sig[j] + 2*((1-l/(lag+1))*gamma(l,t[j], end))
     }
   }
  }
  
  # return list
  return(list(time = t, sig = sig))
}


est.sigma.raw <- function(data, hd, kern, t.index, t.points){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  # mode-handling
  mode = NA
  if(missing(t.index) & missing(t.points)){
    mode = 1
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind = 1:(length(data$time))
  }
  else if(missing(t.index) & !missing(t.points)){
    mode = 2
    t<-t.points
    ind = numeric(length(t))
    for(i in 2:length(t)){
      ind[i] = which.max(data$time[data$time<t[i]])
    }
  }
  else{
    mode = 3
    t<-data$time[t.index]
    ind = t.index
  }
  #t should now be data$time points
  
  tt = length(t)
  n = length(data$time)
  mu = numeric(tt)          # We can only have bandwidth to end amount of calcs
  
  dy = diff(data$Y)
  
  # Optimization removed
  if(mode == 1){
    for(j in 1:tt){
      mu[j] = sqrt((1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*(dy[1:(n-1)])^2))   
    }
  }
  else if(mode == 3){
    for(j in 1:tt){
      mu[j] = sqrt((1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*(dy[1:(n-1)])^2))
    }
  }
  else{ # time points not implemented
    mu[j] = NA
  }
  return(list(time = t, sigRaw = mu))
}

# faster versions
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
    # husk at j bruges i udregning 1 - husk at kør j<-1 hver gang der tjekkes op mod noget i sig[1]!!!!
    
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
    tempData<-est.EveryOtherDiffData(data) 
    data<-NULL #Handles issues with data.frame input
    data$time<-tempData$time #compatability
    dy<-tempData$dy #compatability
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

est.EveryOtherDiffData<-function(data){
  #data as in other estimates functions
  
  diffY<-diff(data$Y)
  everySecondSeq<-seq(1L, length(diffY), by=2L)
  
  dy<-diffY[everySecondSeq]
  
  tempTime <- data$time[seq(1L, length(data$time), by = 2L)] 
  if(length(tempTime)<length(dy)){ #handles uneven vs even number of input
    tempTime<-c(tempTime, 0) #0 unused, just need the right length
  }
  
  return(list(time = tempTime, dy = dy))
}
