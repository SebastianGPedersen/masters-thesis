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


tstar<-function(data, ngroups, trun=floor){
  # data should contain time and test
  out <- apply.partition(data$test, ngroups, data$time, max, trun)
  return(list(start = out$start, end = out$end, tstar = out$res))
}

apply.partition<-function(x, ngroups, time=NA, func = max, trun=floor){
  # x should be vector
  n <- length(x)
  k <- ngroups
  res<-numeric(k)
  
  # if time is not given, do 1:n
  if(missing(time)) time<-1:n
  
  if((n/k)%%1 > 0){
    warning("Imperfect partitioning of data (n/k): ", n/k)     # Lav langt bedre warning!!
  }
    for(j in 1:k){
    #print((trun(n/k)*(j-1)+1):(trun(n/k)*j)) # for help
    res[j] <- func(x[(trun(n/k)*(j-1)+1):(trun(n/k)*j)])  #(floor(n/k)*(j-1)+1):(floor(n/k)*j)
  }
  
  # returns list of the grouped and the raw
  return(list(start = time[(trun(n/k)*((1:k)-1)+1)], end = time[(trun(n/k)*(1:k))], res = res) )
}











# FUTURE

# Calculate Z

#Zcorr = numeric(m)
#for(i in 1:m){
#  fit <- ar(t[i,], order.max = 1) # order.max = 1 corresponds to AR(1)
#  Zcorr[i] <- fit$ar
#}

# Simulate Z from fitted

# Return Z*
