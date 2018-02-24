teststat<-function(data.mu, data.sig, ht, kern){
  # data.mu should include time / mu
  # data.sig should include time / sig
  # kern MUST contain ksq as it is being used
  
    # if whole list is passed - it will find ksq itself
  if(is.list(kern)) ksq = kern$ksq
  if(is.function(kern)) stop("kern argument should be either ksq or list containing ksq")
  
  time <- data.mu$time
  mu<-data.mu$mu
  time2 <- data.sig$time
  sig <- data.sig$sig
  
  # check if lengths makes sense
  if(length(time) != length(mu)){
    stop("length(time) != length(mu)")
  }
  if(length(time2) != length(sig)){
    stop("length(mu) != length(sig)")
  }
  
  # check if time and time2 are NOT identical
  if(!setequal(time, time2)){
    print("time.mu and time.sig differs - attempting to match")
    
    int<-intersect(time,time2)
    #adjust time, mu and sig
    mu<-mu[match(int, time)]
    sig<-sig[match(int, time2)]
    time<-int
  }
  if(length(mu) == length(sig)) print("lengths succesfully matched")
  else stop("Unable to match - clean up your input")
  
  n <- length(time)
  coef <- sqrt(ht/ksq)
  
  # calculates t-stat
  t <- coef*(mu/sig)
  
  # returns list of the grouped and the raw
  return(list(time = time, test = t))
}


# change this
tstar<-function(data, ngroups, trun=floor){
  # data should contain time and test
  # Need vector/list of start and end for each period (or n. of observations in)
  start = data$time[1]
  end = data$time[length(data$time)]
  tstar = max(data$test)
  
  return(list(start = out$start, end = out$end, tstar = out$res))
}

