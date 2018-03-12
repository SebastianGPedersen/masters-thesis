setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")

setting <- sim.setup(Npath = 2, Nstep = 10000)
sims<-sim.heston(setting)
sim<-sim.path(1, sims)

tind<-seq(from=2, to=10000, by=100)

est.mu2 <- function(data, hd, kern, t.index=NA, t.points=NA){
  # data list should include a times column and the Y column (log returns)
  # t.index should be index - we use data$time[t]
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  
  n = length(data$time)
  
  # mode-handling
  mode = NA
  if(is.na(t.index) & is.na(t.points)){
    mode = 1
    t<-data$time[1:(length(data$time))] # if nothing specified - every point in data
  }
  else if(is.na(t.index) & !is.na(t.points)){
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
    if(min(t.index) < 2){
      stop("t.index cannot be less than two (deltaY[1] is not defined)")
    }
    ind<-t.index
  }
  # t should now be data$time points
  
  tt = length(t)
  mu = numeric(tt)
  
  dy = diff(data$Y) # dy without zero - dy[1] = y

  if(mode == 1){
    for(j in 1:(tt-1)){
      mu[j] = (1/hd)*sum(kern((data$time[1:j] - t[j])/hd)*dy[1:j])   
    }
    mu[tt] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[tt])/hd)*dy[1:(n-1)])
    # last point breaks our beautiful method 
  }
  else if(mode == 3){
    for(j in 1:tt){
      mu[j] = (1/hd)*sum(kern((data$time[1:(ind[j])] - t[j])/hd)*dy[1:ind[j]])
    }
    mu[tt] = (1/hd)*sum(kern((data$time[1:(n-1)] - t[tt])/hd)*dy[1:(n-1)])
    
  }
  else{ # time points not implemented
    mu[j] = NA
  }
  
  # return list
  return(list(time = t, mu = mu))
}

est.sigma2 <- function(data, hv, kern, wkern, t.index=NA, lag="auto"){   # we could do lag = "auto"
  # data list should include a times column and the Y column (log returns)
  
  # Handle lag
  if(lag=="auto") lag = 15 #temp
  
  # kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  n = length(data$time)
  # if missing handling:
  if(is.na(t.index)){
    start = lag+1
    t<-data$time[start:n] # if nothing specified - every point in data
    ind<-start:n
  }else{
    t<-data$time[t.index]
    ind<-t.index
  }
  # t should now be data$time points
  
  tt = length(t)
  sig = numeric(tt) 
  
  #dy = cbind(0,t(diff(t(data$Y))))               # diff only does each column seperately / so we transpose to get row wise
  dy = diff(data$Y)
  
  gamma2<-function(l, t){
    # to be used in the loop with j as n            
    {
      # if j+l+1 < n then we should not hit out of range
      #if(ind[j]+l < n){
      #if(1 > 2){
      #  out<- sum(   kern( (data$time[(l+1):(ind[j]+l)] - t)/hv )* 
      #                 dy[(l+1):(ind[j]+l)]*       
      #                 kern( (data$time[1:ind[j]] - t)/hv )*   # TEGN Hvordan Autocov vinduet ser ud
      #                 dy[1:ind[j]]   )
      #  if(is.na(out)) stop("out is NA")
      #  
      #}
      #else{
      #  out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )* 
      #                 dy[(l+1):(n-1)]*       #indices are literally #1 reason for bugs
      #                 kern( (data$time[1:(n-l-1)] - t)/hv )*
      #                 dy[1:(n-l-1)]   )
      #}
    }
    
    out<-ifelse(ind[j]+l < n, sum(   kern( (data$time[(l+1):(ind[j]+l)] - t)/hv )* 
                                  dy[(l+1):(ind[j]+l)]*       
                                  kern( (data$time[1:ind[j]] - t)/hv )*   # TEGN Hvordan Autocov vinduet ser ud
                                  dy[1:ind[j]]   ),
                        sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )* 
                                  dy[(l+1):(n-1)]*       #indices are literally #1 reason for bugs
                                  kern( (data$time[1:(n-l-1)] - t)/hv )*
                                  dy[1:(n-l-1)]   ))
    
    return(out)
  }
  
  gamma<-function(l, t){
    # the standard where end = n
    l <- abs(l)
    out<- sum(   kern( (data$time[(l+1):(n-1)] - t)/hv )* 
                   dy[(l+1):(n-1)]*       #indices are literally #1 reason for bugs
                   kern( (data$time[1:(n-l-1)] - t)/hv )*
                   dy[1:(n-l-1)]   )
    return(out)
  }
  
  for (j in 1:(tt-1)) {
    sig[j] = sum(  (kern(   (data$time[1:ind[j]] - t[j])/hv   )*dy[1:ind[j]])^2  )  # l = 0
    if (lag >=1) {
      for(l in 1:lag){
        sig[j] = sig[j] + 2*(wkern(l/(lag+1))*gamma2(l,t[j])) #gamma2
                                                                   
      }
    }
  }
  sig[tt] = sum(  (kern(   (data$time[1:(n-1)] - t[tt])/hv   )*dy[1:(n-1)])^2  )  # l = 0
  if (lag >=1) {
    for(l in 1:lag){
      sig[tt] = sig[tt] + 2*(wkern(l/(lag+1))*gamma(l, t[tt]))
    }
  }
  
  # return list
  return(list(time = t, sig = sig))
}
tind

(test<-est.mu(data = sim, hd = 0.001, kern = kern.leftexp, t.index = tind)$mu-est.mu2(data = sim, hd = 0.001, kern = kern.leftexp, t.index = tind)$mu)

(lags<-(est.sigma(lag=0, data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen)$sig
       -est.sigma2(lag=0, data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen)$sig))

#(test<-(est.sigma(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen)$sig
#       -est.sigma2(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen)$sig))

(netto<-(est.sigma(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen, t.index = tind)$sig
         -est.sigma2(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen, t.index = tind)$sig))

require(microbenchmark)

microbenchmark(est.sigma(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen, t.index = tind),
               est.sigma2(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen, t.index = tind), times = 10)

require(profvis)
profvis(est.sigma(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen))
profvis(est.sigma2(data = sim, hv = 0.001, kern = kern.leftexp, wkern = kern.parzen))