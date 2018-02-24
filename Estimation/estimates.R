est.mu <- function(data, hd, kern, t.index=NA, t.points=NA){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  # mode-handling
  mode = NA
  if(missing(t.index) & missing(t.points)){
    mode = 1
    t<-1:(length(data$time)) # if nothing specified - every point in data
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

  n = length(t)
  mu = numeric(n)          # We can only have bandwidth to end amount of calcs
  
  #dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  dy = c(0,diff(data$Y))
  # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  
  
  if(mode == 1){
    for(j in 1:n){
      mu[j] = (1/hd)*sum(kern((data$time[1:(j-1)] - t[j])/hd)*dy[2:j])   
    }
  }
  else{
    for(j in 1:n){
      mu[j] = (1/hd)*sum(kern((data$time[1:ind[j]] - t[j])/hd)*dy[2:n])   
    }
  }
  
  #for(j in 1:n){
  #  mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[2:n])   
  #}
  
  # return list
  return(list(time = t, mu = mu))
}

est.sigma <- function(data, hv, lag="auto", kern, wkern, t.index=NA, t.points=NA){   # we could do lag = "auto"
  # data list should include a times column and the Y column (log returns)
  
  # Handle lag
  if(lag=="auto") lag = 10 #temp
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  # if missing handling:
  if(missing(t.index) & missing(t.points)){
    start = lag+1
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
  
  #dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  dy = c(0,diff(data$Y))
  # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  
  gamma<-function(l, t){
    l <- abs(l)
    out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )*dy[(1+l+1):n]*       #indices are literally #1 reason for bugs
                   kern( (data$time[1:(n-l-1)] - t ) * dy[2:(n-l)]   ) )  
    return(out)
  }
  
  L = -lag:lag

  for(j in 1:n){
    for(l in L){                   # Optimize this?
      sig[j] = sig[j] + wkern(l/n)*gamma(l,t[j])
    }
  }
  # return list
  return(list(time = t, sig = sig))
}




