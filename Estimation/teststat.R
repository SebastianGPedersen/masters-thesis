teststat<-function(time, mu, sig, ht, kern, time2=NA){
  # list or individual input?
  # kern MUST contain ksq as it is being used
  
    # if whole list is passed - it will find ksq itself
  if(is.list(kern)) ksq = kern$ksq
  if(is.function(kern)) stop("kern argument should be either ksq or list containing ksq")
  
  # check that time mu and sig has same length
  if(missing(time2)){
    if(length(time) != length(mu)){
      stop("length(time) != length(mu)")
    }
    if(length(mu) != length(sig)){
      stop("length(mu) != length(sig)")
    }
  }
  else{
    # if time2 is given, we will only calculate on intersection (snittet)
     time <- intersect(time,time2)
     
     # goodluck finding indices that corresponds to the new time
     # mu = ~~~~ 
     # sig = ~~~~               #overwrites
  }
  
  n <- length(time)
  coef <- sqrt(ht/ksq)
  
  # calculates t-stat
  t <- coef*(mu/sig)
  
  # returns list of the grouped and the raw
  return(list(time = time, test = t))
}


tstar<-function(data, ngroups, trun=floor){
  # data should contain time and test
  apply.partition(data$test, ngroups, data$time, func=max, trun=floor)
}

apply.partition<-function(x,ngroups, time=NA, func = max, trun=floor){
  # x should be vector
  n <- length(x)
  k <- ngroups
  tstar<-numeric(k)
  
  # if time is not given, do 1:n
  if(missing(time)) time<-1:n
  
  if(n/k%%1 > 0){
    warning("Imperfect partitioning of data (n/k):" +n/k)     # Lav langt bedre warning!!
  }
    for(j in 1:k){
    #print((trun(n/k)*(j-1)+1):(trun(n/k)*j)) # for help
    res[j] <- func(x[(trun(n/k)*(j-1)+1):(trun(n/k)*j)])  #(floor(n/k)*(j-1)+1):(floor(n/k)*j)
  }
  
  # returns list of the grouped and the raw
  return(list(gstart = time[(trun(n/k)*((1:k)-1)+1)], gend = time[(trun(n/k)*(1:k))], res = res) )
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
