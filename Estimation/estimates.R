est.mu <- function(data, hd, kern, t.index=NA, t.points=NA){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # if missing handling:
  if(missing(t.index) & missing(t.points)){
    t<-1:(length(data$time)) # if nothing specified - every point in data
  }
  else if(missing(t.index) & !missing(t.points)){
    t<-t.points
  }
  else{
    t<-data$time[t.index]
  }
  #t should now be data$time points

  n = length(t)
  mu = numeric(n)          # We can only have bandwidth to end amount of calcs
  
  dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  
  for(j in 1:n){
    mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[2:n])   # maybe limit such that it doesnt include                                                                                    
  }                                                                                      # future values that become zero anyway (1:t)
  # return list
  return(list(time = t, mu = mu))
}

est.sigma <- function(data, hv, lag, kern, wkern, t.index=NA, t.points=NA){
  # data list should include a times column and the Y column (log returns)
  
  start = lag+1
  # if missing handling:
  if(missing(t.index) & missing(t.points)){
    t<-start:(length(data$time)) # if nothing specified - every point in data
  }
  else if(missing(t.index) & !missing(t.points)){        # SET UP WARNINGS IF POINTS/INDEX ARE NOT SMART
    t<-t.points
  }
  else{
    t<-data$time[t.index]
  }
  # t should now be data$time points
  
  n = length(t)
  sig = numeric(n)         
  
  dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  
  gamma<-function(l, t){
    l <- abs(l)
    out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )*dy[(1+l+1):n]*       #indices are literally #1 reason for bugs
                   kern( (data$time[1:(n-l-1)] - t ) * dy[2:(n-l)]   ) )  
    return(out)
  }
  
  L = -lag:lag

  for(j in 1:n){
    for(l in L){                   # Optimize this!!!!
      sig[j] = sig[j] + parzenkern(l/n)*gamma(l,t[j])
    }
  }
  # return list
  return(list(time = t, sig = sig))
}




