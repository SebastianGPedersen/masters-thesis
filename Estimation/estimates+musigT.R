estimates <- function(data, hd, hv, lag, kern){
  # data list should include a times column and the Y column (log returns)
  n = length(data$time)
  m = dim(data$Y)[1]                             # should very well be Npaths
  start = lag+1
  mu = matrix(0, nrow = m, ncol = n-start)          # We can only have bandwidth to end amount of calcs
  sig = matrix(0, nrow = m, ncol = n-start)         
  
  dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  
  gamma<-function(l, i, t){
    l <- abs(l)
    out<- sum(   kern( (data$time[(l+1):(n-1)] - data$time[t])/hv )*dy[i,(1+l+1):n]*       #indices are literally #1 reason for bugs
                   kern( (data$time[1:(n-l-1)] - data$time[t] ) * dy[i,2:(n-l)]   ) )  
    # î--- fix this - 100p fejl somewhere
    return(out)
  }
  
  L = -lag:lag
  #mu estimation
  
  #pre calc kerns - very repetitive form - easily save comp time if we use indices in a smart way
  prekern = matrix(0, nrow = n, ncol = n-1)
  for(t in start:n){
    prekern[t,] <- kern((data$time[1:(n-1)] - data$time[t])/hd)
  }
  
  # repeat with the gamma's kernels?
  
  for(i in 1:m)
  {
    #Can we remove the t loop in some way?
    for(t in start:n){
      #mu[i,t-start] = (1/hd)*sum(kern((data$time[1:(n-1)] - data$time[t])/hd)*dy[i,2:n])   # maybe limit such that it doesnt include
      mu[i,t-start] = (1/hd)*sum(prekern[t,]*dy[i,2:n])                                                                       # future values that become zero anyway (1:t)
    }
    
    # sig
    for(t in start:n){
      for(l in L){                   # Optimize this!!!!
        sig[i,t-start] = sig[i, t-start] + parzenkern(l/n)*gamma(l,i,t)
      }
    }
  }
  # return list
  return(list(time = data$time[start:n], mu = mu, sig = sig))
}

teststat<-function(estimates, ht, ksq, ngroups){
  # estimates list should include time, mu, sig
  # possible to generalize tstar handling by inputing function (max)
  
  n <- length(estimates$time)
  m = dim(estimates$mu)[1]
  k <- ngroups
  tstar<-matrix(0, nrow = m, ncol = k)
  
  if(n/k%%1 > 0){
    print("Imperfect partitioning of data (n/k)")     # Lav langt bedre warning!!
  }
  coef <- sqrt(ht/ksq)
  
  
  # calculates t-stat
  t <- coef*(estimates$mu/estimates$sig)
  for(i in 1:m){
    for(j in 1:k){
      print((floor(n/k)*(j-1)+1):(floor(n/k)*j))
      tstar[i,j] <- max(t[i,(floor(n/k)*(j-1)+1):(floor(n/k)*j)])  #(floor(n/k)*(j-1)+1):(floor(n/k)*j)
    }
  }
  
  # Calculate Z
  
  #Zcorr = numeric(m)
  #for(i in 1:m){
  #  fit <- ar(t[i,], order.max = 1) # order.max = 1 corresponds to AR(1)
  #  Zcorr[i] <- fit$ar
  #}
  
  # Simulate Z from fitted
  
  # Return Z*
  
  
  # returns list of the grouped and the raw
  return(list(grouped = list(gstart = (floor(n/k)*((1:k)-1)+1), gend = (floor(n/k)*(1:k)), tstar = tstar) ), raw = list(time = estimates$time, test = t) )
}