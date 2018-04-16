
# BANDWIDTH SHOULD BE TRANSLATED FROM SECONDS TO YEARS ( BW / Seconds per year)
# Consider multiplying time such that our unit is in seconds and not years...!

est.mu <- function(data, hd, kern = kern.leftexp, t.index, t.points){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # t.points does NOT work yet!
  
  #if(is.null(data$dy)){
  #  dy <- diff(data$Y)
  #} else {
  #  dy <- data$dy
  #}
  dy <- data$Y
  
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

est.sigma <- function(data, hv, kern = kern.leftexp, wkern = kern.parzen, t.index, lag="auto"){   # we could do lag = "auto"
  # data list should include a times column and the Y column (log returns)
  
  lagset <- lag
  
  #if(is.null(data$dy)){
  #  dy <- diff(data$Y)
  #} else {
  #  dy <- data$dy
  #}
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
    # more elegant solution to nmu should be implemented
    nmu <- numeric(tt)
    for(j in 1:tt){
      nmu[j] <- sum(kern((data$time[1:(n-1)] - t[j])/hd))
    }
  }
  
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
    if(lagset =="auto")
    {
      lag <- laglength(dy, nmu[j])       # NOT EXACTLY LIKE KIM'S / BUT SHOULD NOT MAKE A HUGE ENOUGH IMPACT
      #print(lag)
    }
   sig[j] <- sum(  (kern(   (data$time[1:(end-1)] - t[j])/hv)*dy[1:(end-1)]  )^2  )  # l = 0
   
   #sig[j] = sum (  dy[2:end]^2  )
   if (lag >=1) {
     for(l in 1:lag){
       #sig[j] = sig[j] + 2*(wkern(l/n)*gamma(l,t[j], end))
       sig[j] <- sig[j] + 2*(wkern(l/(lag))*gamma(l,t[j], end))
       #sig[j] = sig[j] + 2*((1-l/(lag+1))*gamma(l,t[j], end))
     }
   }
  }
  sig <- sig/hv # NEW FROM KIM'S
  
  # return list
  return(list(time = t, sig = sig))
}

est.sigma.raw <- function(data, hv, kern, t.index, t.points){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  # mode-handling
  mode <- NA
  if(missing(t.index) & missing(t.points)){
    mode <- 1
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
    ind <- 1:(length(data$time))
  }
  else if(missing(t.index) & !missing(t.points)){
    mode <- 2
    t<-t.points
    ind <- numeric(length(t))
    for(i in 2:length(t)){
      ind[i] <- which.max(data$time[data$time<t[i]])
    }
  }
  else{
    mode <- 3
    t<-data$time[t.index]
    ind <- t.index
  }
  #t should now be data$time points
  
  tt <- length(t)
  n <- length(data$time)
  sig <- numeric(tt)          # We can only have bandwidth to end amount of calcs
  
  dy <- data$Y
  
  # Optimization removed
  if(mode == 1){
    for(j in 1:tt){
      sig[j] <- sqrt((1/hv)*sum(kern((data$time[1:(n-1)] - t[j])/hv)*(dy[1:(n-1)])^2))   
    }
  }
  else if(mode == 3){
    for(j in 1:tt){
      sig[j] <- sqrt((1/hv)*sum(kern((data$time[1:(n-1)] - t[j])/hv)*(dy[1:(n-1)])^2))
    }
  }
  else{ # time points not implemented
    sig[j] <- NA
  }
  return(list(time = t, sig = sig))
}

# faster versions
est.mu.next <-function(data, prevmu, hd, t.index){
  # data should contain time | logreg - prevmu should contain time | mu
  # t.index indicates where we wish to calculate the estimator
  
  # --- debug ---
  # data<-list(time = 1:10/10, Y = 1:10)
  # hd = 0.1
  # prevmu <- est.mu(data, hd, t.index = 1:3, originalEstimator = T)
  # t.index<-c(3,5,9)
  
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
  sig[1] <- startsig$sig*rescale[1] + (1/hv)*sum(  (kern(   (data$time[start[1]:end[1]] - t[1])/hv   )*dy[start[1]:end[1]])^2 )  # l = 0
  if (lag >=1) {
    for(l in 1:lag){
      sig[1] <- sig[1] + 2*(1/hv)*(wkern(l/(lag))*gamma(l) )
    }
  }
  for(j in 2:tt){
    sig[j] <- sig[j-1]*rescale[j] + (1/hv)*sum(  (kern(   (data$time[start[j]:end[j]] - t[j])/hv   )*dy[start[j]:end[j]])^2  )  # l = 0
    if (lag >=1) {
      for(l in 1:lag){
        sig[j] <- sig[j] + 2*(1/hv)*(wkern(l/(lag))*gamma(l) )
      }
    }
  }
  
  # --- debug ---
  #{
    # husk at j bruges i udregning 1 - husk at kør j<-1 hver gang der tjekkes op mod noget i sig[1]!!!!
    
    #(1/hv)*(sum( (kern( (data$time[1:(n-1)] - t[1])/hv)*dy[1:(n-1)])^2 )+2*wkern(1)*gamma2(1, t[1]))
    #startsig$sig<-(1/hv)*(sum( (kern( (data$time[1:(n-1)] - data$time[4])/hv)*dy[1:(n-1)])^2 )+2*wkern(1)*gamma2(1, data$time[4]))
    
    #gamma2<-function(l, t){
    #  out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )*dy[(l+1):(n-1)]*       
    #                 kern( (data$time[1:(n-l-1)] - t)/hv )*dy[1:(n-l-1)]   )
    #  return(out)
    #}
  #}
  
  #ADD BLOCK TO EXISTING (wow much blockchainy) (if it existed)
  t <- c(prevsig$time, t)
  sig <- c(prevsig$sig, sig) #we should remember to divide with hv
  
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
  
  dy = diff(data$Y)

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

# LAG LENGTH #
# find Newey-West (1994) -- automatic lag selection with parzen kernel

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

  #if(kern == kern.leftexp)
  #{
    # mu.next only works for left exp kern
  out<-est.mu.next(data = data, hd = hd, t.index = t.index)
  #}
  #else{
  #  out <- est.mu(data = data, hd = hd, kern = kern, t.index = t.index)
  #}
  n <- length(data$Y)
  psi1 = 0.25
  coef <- n/(n-kn)*1/(psi1*kn)
  
  return(list(time = out$time, mu = out$mu*coef))
}

est.sigma.new <- function(data, hv, kern = kern.leftexp, t.index, kn, arfun, theta){
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
  
  out <- est.sigma.raw(data = data, hv = hv, kern = kern, t.index = t.index)
  n <- length(data$Y)
  psi2 = 0.25
  coef <- n/(n-kn)*1/(psi2*kn)
  
  #prep args
  args <- list(data = data, hv = hv, kern = kern, t.index = t.index)
  
  ar<-arfun(args, theta)
  
  return(list(time = out$time, sig = out$sig^2*coef-ar))
}

est.ar.iid <- function(args, theta){
  # calculates omega^2
  # data needs to include un-preaveraged
  # args is a list that contains all the info that the ar needs
  data <- args$data
  hv <- args$hv
  kern <- args$kern
  t.index <- args$t.index
  
  psi1 <- 0.25
  psi2 <- 1/12
  
  dy <- diff(data$raw)
  
  t<-data$time[t.index]
  ind = t.index

  tt = length(t)
  n = length(data$time)
  omega = numeric(tt)          # We can only have bandwidth to end amount of calcs
  
  for(j in 1:tt){
    omega[j] = (1/hv)*sum(kern((data$time[2:(n-1)] - t[j])/hd)*dy[2:(n-1)]*dy[1:(n-2)])   
  }
  return(psi1/(psi2*theta^2)*omega)
}

