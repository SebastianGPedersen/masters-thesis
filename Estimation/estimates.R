source("kernels/kernels.R")
# BANDWIDTH SHOULD BE TRANSLATED FROM SECONDS TO YEARS ( BW / Seconds per year)
# Consider multiplying time such that our unit is in seconds and not years...!

est.mu <- function(data, hd, kern = kern.leftexp, t.index){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  n = length(data$time)
  
  if(missing(t.index)){
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind<-1:length(data$time)
  }
  else{
    t<-data$time[t.index]
    ind<-t.index
  }
  # t should now be data$time points
  
  dy = data$Y
  
  #update n accordingly
  n = length(data$time)
  tt = length(t)
  mu = numeric(tt)
  
  for(j in 1:(tt)){
    mu[j] = (1/hd)*sum(kern((data$time[1:(ind[j])] - t[j])/hd)*dy[1:ind[j]])
  }
  
  # return list
  return(list(time = t, mu = mu))
}

est.sigma <- function(data, hv, t.index, kern = kern.leftexp, wkern=kern.parzen, lag="auto"){
  # data list should include a times column and the Y column (log returns)
  
  lagset <- lag
  dy <- data$Y
  
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
  
  # Handle lag
  if(lagset=="auto"){
    nmu <- numeric(tt)
    for(j in 1:tt){
      nmu[j] <- sum(kern((data$time[1:(n-1)] - t[j])/hv))
    }
  }
  
  gamma<-function(l, t){
    # end is the highest needed index
    out<- sum(   kern( (data$time[(l+1):ind[j]] - t)/hv )* 
                   dy[(1+l):(ind[j])]*       
                   kern( (data$time[1:(ind[j]-l)] - t)/hv )*
                   dy[1:(ind[j]-l)]   )
    return(out)
  }
  
  end = n
  for (j in 1:tt) {
    if(lagset =="auto"){
      lag <- laglength(dy, nmu[j])
    }
    
    sig[j] <- sum(  (kern(   (data$time[1:ind[j]] - t[j])/hv)*dy[1:ind[j]]  )^2  )  # l = 0
    if (lag >=1) {
      for(l in 1:lag){
        sig[j] <- sig[j] + 2*(wkern((l-1)/(lag))*gamma(l,t[j]))
      }
    }
  }
  sig <- sig/hv
  
  # return list
  return(list(time = t, sig = sig))
}

est.sigma.raw <- function(data, hv, kern = kern.leftexp, t.index){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  # Calculates sigma^2
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  # mode-handling
  if(missing(t.index)){
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind <- 1:(length(data$time))
  }
  else{
    t<-data$time[t.index]
    ind <- t.index
  }
  #t should now be data$time points
  
  tt <- length(t)
  n <- length(data$time)
  sig <- numeric(tt)
  
  dy <- data$Y
  
  # Optimization removed
  for(j in 1:tt){
    sig[j] <- (1/hv)*sum(kern((data$time[1:ind[j]] - t[j])/hv)*(dy[1:ind[j]])^2)
  }
  return(list(time = t, sig = sig))
}

# faster versions
est.mu.next <-function(data, prevmu, hd, t.index){
  # data should contain time | logreg - prevmu should contain time | mu
  # t.index indicates where we wish to calculate the estimator
  
  # --- debug ---
  {
    # data<-list(time = 1:10/10, Y = 1:10)
    # hd = 0.1
    # prevmu <- est.mu(data, hd, t.index = 1:3, originalEstimator = T)
    # t.index<-c(3,5,9)
  }
  
  dy <- data$Y
  
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
  
  if(anyNA(dy[t.index])) warning("t.index cannot be higher than t_n-1 (this will rightfully give you NAs)")
  
  n <- length(data$time)
  
  # init
  tt = length(t)
  mu <- numeric(tt)
  # CALCULATE NEXT 'BLOCK'
  mu[1] <- startmu$mu*rescale[1] + 1/hd * sum( kern( (data$time[start[1]:end[1]] - t[1])/hd)*dy[start[1]:end[1]] )
  if(tt > 1){
    for(j in 2:(tt)){
      mu[j] <- mu[j-1]*rescale[j] + 1/hd * sum( kern(  (data$time[start[j]:end[j]] - t[j])/hd)*dy[start[j]:end[j]] )
    }
  }
  
  # debug checker
  #1/hd * sum( kern( (data$time[1:(n-1)] - data$time[5])/hd)*dy[1:(n-1)] )
  
  #ADD BLOCK TO EXISTING (wow much blockchainy) (if it existed)
  #if(min(prevmu$time) >= 0){
  #  t <- c(prevmu$time, t)
  #  mu <- c(prevmu$mu, mu)
  #} 
  
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
  
  dy <- data$Y
  
  # --- debug ---
  {
    #data<-list(time = 1:100/100, Y = 1:100)
    #hv = 0.1
    #prevsig <- est.sigma(data, hv=hv, t.index = 2:4, lag = 1)
    #t.index<-c(6,9)
  }
  
  if(missing(prevsig)){
    # initial handling is pretty weird because of lag length
    prevsig <- est.sigma(data, hv=hv, t.index = t.index[1],lag = lag) # change to minus
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
  sig[1] <- startsig$sig*rescale[1] + (1/hv)*sum(  (kern(   (data$time[start[1]:end[1]] - t[1])/hv   )*dy[start[1]:end[1]])^2 )  # l = 0
  if (lag >=1) {
    for(l in 1:lag){
      sig[1] <- sig[1] + 2*(1/hv)*(wkern(l/(lag))*gamma(l) )
    }
  }
  if(tt >= 2){
    for(j in 2:tt){
      sig[j] <- sig[j-1]*rescale[j] + (1/hv)*sum(  (kern(   (data$time[start[j]:end[j]] - t[j])/hv   )*dy[start[j]:end[j]])^2  )  # l = 0
      if (lag >=1) {
        for(l in 1:lag){
          sig[j] <- sig[j] + 2*(1/hv)*(wkern(l/(lag))*gamma(l) )
        }
      }
    }
  }
  
  # --- debug ---
  {
    # husk at j bruges i udregning 1 - husk at kør j<-1 hver gang der tjekkes op mod noget i sig[1]!!!!
    
    #(1/hv)*(sum( (kern( (data$time[1:(n-1)] - t[1])/hv)*dy[1:(n-1)])^2 )+2*wkern(1)*gamma2(1, t[1]))
    #startsig$sig<-(1/hv)*(sum( (kern( (data$time[1:(n-1)] - data$time[4])/hv)*dy[1:(n-1)])^2 )+2*wkern(1)*gamma2(1, data$time[4]))
    
    #gamma2<-function(l, t){
    #  out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )*dy[(l+1):(n-1)]*       
    #                 kern( (data$time[1:(n-l-1)] - t)/hv )*dy[1:(n-l-1)]   )
    #  return(out)
    #}
  }
  
  #ADD BLOCK TO EXISTING (wow much blockchainy) (if it existed)
  t <- c(prevsig$time, t)
  sig <- c(prevsig$sig, sig) #we should remember to divide with hv
  
  # return list
  return(list(time = t, sig = sig))
}

est.mu2 <- function(data, hd, kern = kern.leftexp, t.index, t.points){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  n = length(data$time)
  
  # mode-handling
  mode = 0
  if(missing(t.index)){
    mode = 1
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind<-1:length(data$time)
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
  
  dy = data$Y

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

laglength = function(dx, nmu){
  # dx are NOT preavr
  # only dx from backwards in time should be considered
  # Indices may be wrong due to zero index (or is it? Matlab is index 1...)
  c <- 2.6614
  q <- 2
  
  n <- round(4*(nmu/100)^(4/25))
  root <- 1/(2*q+1)
  
  # pre-estimate autocovariance of noise increment
  E_deij <- numeric(n+1)
  for(i in 1:(n+1)){
    E_deij[i] <- mean(dx[(i+1):length(dx)]*dx[1:(length(dx)-i)])
  }
  
  v0 <-  -sum(   (1:(n+1))*E_deij   )
  
  ac<- numeric(n)
  for(i in 1:n){
    ac[i] <- -sum(   (1:(n+1-i))*E_deij[(i+1):(n+1)]   )
  }
  
  # determine lag length
  
  s0 <- v0 + 2 * sum(ac)
  sq <- 2*sum(   (1:n)^(q)*ac   ) # sq = 2*sum((1:n).^q.*ac); - regnehieraki
  
  gamma <- c*(   ((sq/s0)^q)^root   )
  return(round(gamma*nmu^root))
}

est.mu.new <- function(data, hd, kern = kern.leftexp, t.index, kn){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  # This function is more like a wrapper around original mu but with new coefficient
  
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")

  if(identical(kern,kern.leftexp$kern) & length(t.index)>1)
  {
    # mu.next only works for left exp kern
    out<-est.mu.next(data = data, hd = hd, t.index = t.index)
  }
  else{
    out <- est.mu(data = data, hd = hd, kern = kern, t.index = t.index)
  }
  psi1 = 0.25
  coef <- 1/(psi1*kn)
  
  return(list(time = out$time, mu = out$mu*coef))
}

est.sigma.new <- function(data, hv, kern = kern.leftexp, t.index, kn, noisefun, theta){
  # data list should include a time column , Y column (log returns) and non-preavr obs (raw)
  # t.index should be index - we use data$time[t]
  # Returns sigma^2
  
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  {
    #if(kern = kern.leftexp & lag!= auto)
    #{
    #  # mu.next only works for left exp kern
    #  out<-est.sig.next(data = data, hd = hv, t.index = t.index)
    #}
  }
  
  if(identical(kern, kern.leftexp$kern) & length(t.index)>1){
    #print("using .next")
    out <- est.sigma.raw.next(data = data, hv = hv, t.index = t.index)
  }
  else{
    out <- est.sigma.raw(data = data, hv = hv, kern = kern, t.index = t.index)
  }
  
  n <- length(data$Y)
  psi2 <- 1/12
  coef <- 1/(psi2*kn)
  
  #prep args
  args <- list(data = data, hv = hv, kern = kern, t.index = t.index)
  #args <- list(data = data, hv = 300000, kern = kern, t.index = t.index)
  
  noise<-noisefun(args, theta)
  
  return(list(time = out$time, sig = out$sig*coef-noise, noise = noise, sig2 = out$sig*coef))
}

est.noise.iid <- function(args, theta){
  # calculates omega^2
  # data needs to include un-preaveraged
  # args is a list that contains all the info that the ar needs
  data <- args$data
  hv <- args$hv
  kern <- args$kern
  t.index <- args$t.index
  
  psi2 <- 1/12 #int(g(u)^2)
  psi3 <- 1    #int(g'(u)^2)
  
  dy <- diff(data$raw)
  
  t<-data$time[t.index]
  ind <- t.index

  dt <- data$time[2] - data$time[1]
  
  tt <- length(t.index)
  n <- length(data$time)
  omega <- numeric(tt)          # We can only have bandwidth to end amount of calcs
  
  for(j in 1:tt){
    omega[j] <- -(dt/hv)*sum(kern((data$time[2:ind[j]] - t[j])/hv)*dy[2:ind[j]]*dy[1:ind[j]-1]) # remove optimization if this goes badly...
  }
  return(psi3/(psi2*theta^2)*omega)
  #return(omega)
}

est.noise.iid.next <- function(args, theta){
  # calculates omega^2
  # data needs to include un-preaveraged
  # args is a list that contains all the info that the ar needs
  data <- args$data
  hv <- args$hv
  kern <- kern.leftexp$kern
  t.index <- args$t.index
  
  #debug
  {
    #data<-list(time = 1:100/100, Y = 1:100, raw = 1:100)
    #hv = 0.1
    #t.index<-c(6,10,15)
    #kern <- kern.leftexp$kern
    #theta <- 1
  }
  
  dy <- diff(data$raw)
  
  # calculate previous (takes list as argument)
  prevarg <- list(data = data, hv = hv, kern = kern, t.index = t.index[1])
  prev <- est.noise.iid(prevarg, theta)
  startomega<-list(time = data$time[t.index[1]], omega = prev[length(prev)])
  
  if(data$time[t.index][1] < startomega$time) stop("t.index should be higher than previous sig times")
  if(length(t.index) == 1) return(startomega)
  if(data$time[t.index][1] == startomega$time) t.index <- t.index[2:length(t.index)]
  
  # time and index handling
  t<-data$time[t.index]
  end <- t.index
  if(length(t.index) > 1){
    start <- c( min(which(data$time > startomega$time)), end[1:(length(end)-1)]+1)
  }
  else{
    start <- min(which(data$time > startomega$time))
  }
  
  # scaling
  dt1<-t[1]-startomega$time
  rescale <- exp(-c(dt1, diff(t))/hv)
  
  n <- length(data$time)
  
  # THIS SHOULD BE FIXED IN NON-EQUIDISTANT
  dt<-data$time[2]-data$time[1]
  
  psi2 <- 1/12
  psi3 <- 1
  coef <- psi3/(psi2*theta^2)
  # init
  tt = length(t)
  omega <- numeric(tt)
  # CALCULATE NEXT 'BLOCK'
  omega[1] <- startomega$omega*rescale[1] - coef*dt/hv * sum( kern( (data$time[start[1]:end[1]] - t[1])/hv)*dy[start[1]:end[1]]*dy[(start[1]-1):(end[1]-1)] )
  if(tt > 1){
    for(j in 2:(tt)){
      omega[j] <- omega[j-1]*rescale[j] - coef*dt/hv * sum( kern(  (data$time[start[j]:end[j]] - t[j])/hv)*dy[start[j]:end[j]]*dy[(start[j]-1):(end[j]-1)] )
    }
  }
  omega <- c(prev, omega)
  #debug
  {
    #coef*(1/hv)*sum(kern((data$time[2:(n-1)] - t[j])/hv)*dy[2:(n-1)]*dy[1:(n-2)])
  }

  return(omega)
}

est.sigma.raw.next <- function(data, prevsig, hv, t.index){   #
  # data list should include a time column , Y column (log returns) and non-preavr obs (raw)
  # t.index should be index - we use data$time[t]
  # Returns sigma^2
  
  kern <- kern.leftexp$kern
  
  dy <- data$Y
  
  if(missing(prevsig)){
    append = T
    # initial handling is pretty weird because of lag length
    prevsig <- est.sigma.raw(data = data, hv=hv, t.index = t.index[1]) # change to minus
    if(length(t.index) == 1){
      return(prevsig)
    }
  } 
  startsig<-list(time = prevsig$time[length(prevsig$time)], sig = prevsig$sig[length(prevsig$sig)])
  
  if(data$time[t.index][1] < startsig$time) stop("t.index should be higher than previous sig times")
  if(data$time[t.index][1] == startsig$time) t.index <- t.index[2:length(t.index)]
  
  t<-data$time[t.index]
  end <- t.index
  start <- c( min(which(data$time > startsig$time)), end[1:(length(end)-1)]+1)
  # scaling
  dt1<-t[1]-startsig$time
  dt<-c(dt1, diff(t))
  rescale <- exp(-dt/hv)
  
  if(anyNA(dy[t.index])) warning("t.index cannot be higher than t_n-1 (this will rightfully give you NAs)")
  
  n <- length(data$time)
  
  # init
  tt = length(t)
  sig <- numeric(tt)
  
  # CALCULATE NEXT 'BLOCK'
  j<-1
  sig[1] <- startsig$sig*rescale[1] + (1/hv)*sum(  kern(   (data$time[start[1]:end[1]] - t[1])/hv   )*(dy[start[1]:end[1]])^2 )
  if(tt > 1){
    for(j in 2:tt){
      sig[j] <- sig[j-1]*rescale[j] + (1/hv)*sum(  kern(   (data$time[start[j]:end[j]] - t[j])/hv   )*(dy[start[j]:end[j]])^2  )
    }
  }
  
  #ADD BLOCK TO EXISTING (wow much blockchainy) (if it existed) (fuck - we dont want that)
  if(append){
    t <- c(prevsig$time, t)
    sig <- c(prevsig$sig, sig)
  }
  
  # return list
  return(list(time = t, sig = sig))
}

est.mu.mat <- function(data, hd, t.index, kern = kern.leftexp$kern){
  # data contains TIME and a log-return MATRIX
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  n = length(data$time)
  
  if(missing(t.index)){
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind<-1:length(data$time)
  }
  else{
    t<-data$time[t.index]
    ind<-t.index
  }
  # t should now be data$time points
  
  dy <- data$Y #should now be matrixxxzz
  
  #update n accordingly
  n <- length(data$time)
  path <- dim(data$Y)[1]
  if(is.null(path)) stop("Do not use .mat if input is not matrix!")
  tt <- length(t)
  mu <- matrix(NA, nrow = path, ncol = tt)
  
  # calc kerns
  kerns <- matrix(0, nrow = tt, ncol = n)
  for(j in 1:tt){
    kerns[j,1:ind[j]] <- kern(   (data$time[1:(ind[j])] - t[j])/hd)
  }
  
  for(j in 1:(tt)){
    kern_now <- kerns[j,1:ind[j]]
    for(i in 1:path){
      mu[i,j] <- (1/hd)*sum(kern_now*dy[i, 1:ind[j]])
    }
  }
  
  # return list
  return(list(time = t, mu = mu))
}

est.sigma.mat <- function(data, hv, t.index, kern = kern.leftexp, wkern=kern.parzen, lag="auto"){
  # data list should include a times column and the Y column (log returns)
  
  lagset <- lag
  dy <- data$Y
  
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
  path <- dim(data$Y)[1]
  if(is.null(path)) stop("Do not use .mat if input is not matrix!")
  sig = matrix(NA, nrow = path, ncol = tt)
  # Handle lag
  if(lagset=="auto"){
    nmu <- numeric(tt)
    for(j in 1:tt){
      nmu[j] <- sum(kern((data$time[1:(n-1)] - t[j])/hv))
    }
  }
  
  # calc kerns
  kerns <- matrix(NA, nrow = tt, ncol = n) #NA in kerns acts like trip-wires. If wrong index is used, NA error is thrown. perf!
  for(j in 1:tt){
    kerns[j,1:ind[j]] <- kern(   (data$time[1:(ind[j])] - t[j])/hv)
  }
  
  gamma<-function(l){
    # end is the highest needed index
    out<- sum(   kerns[j,(1+l):ind[j]]  *dy[i,(1+l):(ind[j])]*
                 kerns[j,1:(ind[j]-l)]  *dy[i,1:(ind[j]-l)]   )
    return(out)
  }
  
  
  end = n
  for (j in 1:tt) {
    for(i in 1:path){
      if(lagset =="auto"){
        lag <- laglength(dy, nmu[j])
      }
      sig[i,j] <- sum(  (kerns[j, 1:ind[j] ]*dy[i, 1:ind[j]]  )^2  )  # l = 0
      if (lag >=1) {
        for(l in 1:lag){
          sig[i,j] <- sig[i,j] + 2*(wkern(l/(lag))*gamma(l))
        }
      }
    }
  }
  sig <- sig/hv
  
  # return list
  return(list(time = t, sig = sig))
}

est.mu.mat.next <-function(data, hd, t.index){
  # data should contain time | logreg
  # t.index indicates where we wish to calculate the estimator
  
  dy <- data$Y # <- is a ~matrix~
  kern <- kern.leftexp$kern

  t<-data$time[t.index]
  end <- t.index
  
  if(length(t.index) > 1){
    start <- c( which.min(data$time), end[1:(length(end)-1)]+1)
  }
  else{
    start <- which.min(data$time)
  }
  
  # scaling
  dt<-c(0,diff(t))
  rescale <- exp(-dt/hd)
  
  if(anyNA(dy)) warning("NA in dy[, t.index]")
  
  n <- length(data$time)
  path <- dim(data$Y)[1]
  if(is.null(path)) stop("Do not use .mat if input is not matrix!")
  tt <- length(t)
  
  mu <- matrix(NA, nrow = path, ncol = tt)
  
  kerns <- kern( (data$time[start[1]:end[1]] - t[1])/hd)
  # CALCULATE NEXT 'BLOCK'
  for(i in 1:path){
    mu[i,1] <- 1/hd * sum( kerns*dy[i, start[1]:end[1]] ) # mite need loop
  }
  if(tt > 1){
    for(j in 2:(tt)){
      kerns <- kern(   (data$time[start[j]:end[j]] - t[j])/hd   )
      for(i in 1:path){
        mu[i,j] <- mu[i,j-1]*rescale[j] + 1/hd * sum( kerns*dy[i,start[j]:end[j]] )
      }
    }
  }
  
  return(list(time = t, mu = mu))
}

est.sigma.mat.next <- function(data, hv, t.index, wkern=kern.parzen, lag="auto"){   #
  # data list should include a times column and the Y column (log returns)
  
  # Handle lag
  if(lag=="auto") lag = 15 #temp
  if(t.index[1] < lag) stop("t.index can and should NOT be lower than lag length!")
  
  kern <- kern.leftexp$kern
  wkern = kern.parzen$kern
  
  dy <- data$Y
  
  if(length(t.index) == 1){
    startsig <- est.sigma.mat(data, hv=hv, t.index = t.index[1], lag = lag)
    return(startsig)
  }
   
  t<-data$time[t.index]
  end <- t.index
  start <- c( which.min(data$time), end[1:(length(end)-1)]+1)
  
  if(anyNA(dy[,t.index])) warning("NA in dy[,t.index]")
  # scaling
  dt<-c(0, diff(t))
  rescale <- exp(-2*dt/hv)
  
  
  n = length(data$time)
  tt = length(t)
  path <- dim(data$Y)[1]
  if(is.null(path)) stop("Do not use .mat if input is not matrix!")
  if(path)
  sig = matrix(NA, nrow = path, ncol = tt)
  
  gamma<-function(l){
    out<- sum(    kerns[(start[j]):end[j]]    *dy[i,(start[j]):end[j]]*
                  kerns[(start[j]-l):(end[j]-l)]  *dy[i,(start[j]-l):(end[j]-l)]   ) # <-- gamma is wrong!
    return(out)
  }
  {
    #gamma<-function(l){ # sig.next
    #  out<-sum(   kern( (data$time[start[j]:end[j]] - t[j])/hv )    *dy[start[j]:end[j]]*       
    #                kern( (data$time[(start[j]-l):(end[j]-l)] - t[j])/hv )*dy[(start[j]-l):(end[j]-l)]   )
    #  return(out)
    #}
    
    #gamma<-function(l){ # sig.mat
    # out<- sum(   kerns[j,(1+l):ind[j]]  *dy[i,(1+l):(ind[j])]*
    #                kerns[j,1:(ind[j]-l)]  *dy[i,1:(ind[j]-l)]   )
    # return(out)
    #}
  }
  
  # CALCULATE FIRST
  #j<-1 # for gamma(l) which takes j from fct enviroment
  #kerns <- kern(   (data$time[1:(n-1)] - t[1])/hv)
  #for(i in 1:path){
  #  sig[i,1] <- (1/hv)*sum(  (kerns[start[j]:end[j]] *dy[i, start[1]:end[1]])^2 )  # MAKE GAMMAFIRST? ELLER sigma.mat som første??
  #  if (lag >=1) {
  #    for(l in 1:lag){
  #      sig[i,1] <- sig[i,1] + 2*(1/hv)*(wkern(l/(lag))*gamma(l) ) #loop-sum
  #    }
  #  }
  #}
  
  # calculate first sigma
  sig[,1] <- est.sigma.mat(data, hv, t.index[1], lag = lag)$sig
  
  # CALCULATE OTHERS
  if(tt >= 2){
    for(j in 2:tt){
      kerns <- kern(   (data$time[1:(n-1)] - t[j])/hv)
      for(i in 1:path){
        sig[i,j] <- sig[i,j-1]*rescale[j] + (1/hv)*sum(  ( kerns[start[j]:end[j]]*dy[i,start[j]:end[j]])^2  )  # l = 0
        if (lag >=1) { # move if outside loop!
          for(l in 1:lag){
            sig[i,j] <- sig[i,j] + 2*(1/hv)*(wkern(l/(lag))*gamma(l) ) #loop-sum
          }
        }
      }
    }
  }

  # return list
  return(list(time = t, sig = sig))
}
