
# BANDWIDTH SHOULD BE TRANSLATED FROM SECONDS TO YEARS ( BW / Seconds per year)
# Consider multiplying time such that our unit is in seconds and not years...!

est.mu <- function(data, hd, kern, t.index=NA, t.points=NA){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # t.points might not work correctly....
  
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
  
  #dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  dy = c(0,diff(data$Y))
  # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  
  # Optimization removed
  if(mode == 1){
    for(j in 1:tt){
      mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[2:n])   
    }
  }
  else if(mode == 3){
    for(j in 1:tt){
      #mu[j] = (1/hd)*sum(kern((data$time[1:(ind[j]-1)] - t[j])/hd)*dy[ 2:ind[j] ]) # mu[1] = 0...?
      mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[2:n])
    }
  }
  else{ # time points not implemented
    mu[j] = NA
  }
  
  #for(j in 1:n){
  #  mu[j] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[j])/hd)*dy[2:n])   
  #}
  
  # return list
  return(list(time = t, mu = mu))
}

est.sigma <- function(data, hv, kern, wkern, t.index=NA, lag="auto"){   # we could do lag = "auto"
  # data list should include a times column and the Y column (log returns)
  
  # Handle lag
  if(lag=="auto") lag = 15 #temp
  
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
  }
  else{
    t<-data$time[t.index]
    ind<-t.index
  }
  # t should now be data$time points
  
  n = length(data$time)
  tt = length(t)
  sig = numeric(tt) 
  
  #dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  dy = c(0,diff(data$Y))
  # we put in a column of zeros to fit our sizes - dy[i,1] should NEVER be used!
  
  gamma<-function(l, t, end){
    # end is the highest needed index
    l <- abs(l)
    out<- sum(   kern( (data$time[(l+1):(end-1)] - t)/hv )* dy[(1+l+1):end]*       #indices are literally #1 reason for bugs
                 kern( (data$time[1:(end-l-1)] - t)/hv )* dy[2:(end-l)]   )  
    
    # l <- abs(l)
    # out<- sum(    dy[(1+l+1):end]*       #indices are literally #1 reason for bugs
    #              dy[2:(end-l)]   )  
    return(out)
  }
  
  end = n
  for (j in 1:tt) {
   sig[j] = sum(  (kern(   (data$time[1:(end-1)] - t[j])/hv   )*dy[2:end])^2  )  # l = 0
   
   #sig[j] = sum (  dy[2:end]^2  )
   for(l in 1:lag){
     #sig[j] = sig[j] + 2*(wkern(l/n)*gamma(l,t[j], end))
     sig[j] = sig[j] + 2*(wkern(l/(lag+1))*gamma(l,t[j], end))
   }
  }
  
  # return list
  return(list(time = t, sig = sig, deep = deep, parz = parz))
}


