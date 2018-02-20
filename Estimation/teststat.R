teststat<-function(estimates, ht, ngroups){
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
  
  # returns list of the grouped and the raw
  return(list(grouped = list(gstart = (floor(n/k)*((1:k)-1)+1), gend = (floor(n/k)*(1:k)), tstar = tstar) ), raw = list(time = estimates$time, test = t) )
}

teststat<-function(estimates, ht, ngroups){
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
  
  # returns list of the grouped and the raw
  return(list(grouped = list(gstart = (floor(n/k)*((1:k)-1)+1), gend = (floor(n/k)*(1:k)), tstar = tstar) ), raw = list(time = estimates$time, test = t) )
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
