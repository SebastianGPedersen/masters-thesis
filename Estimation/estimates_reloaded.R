#Install required packages if not already installed
packages <- c('foreach')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load packages
invisible(sapply(packages, require, character.only = TRUE))

#Load sources
setwd(Sys.getenv("masters-thesis"))
source("kernels/kernels.R")



######### FUNCTIONS ######### 

### SIGMA ###
#With parallel 
est.sigma.mat.2.0 <- function(data, hv, kern = kern.leftexp, wkern=kern.parzen, lag=10){
  #data <- path
  #hv <- 5 / (60*24*7*52)
  #lag <- 10
  
  #p0 <- Sys.time()
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hv)
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  lags <- lag
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  #temp_func
  sigma_func <- function(path) {
    
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy
    
    #zero lag
    sum_terms <- products^2
    cum_sums <- cumsum(sum_terms)
    gamma_ls <- cum_sums
    sigmas <- gamma_ls[(lags+1):length(gamma_ls)] #It has to fit length-wise with all lags
    
    for (lag in 1:lags) {
      #lag <- 1
      sum_terms <- products[1:(n-lag)]*products[(lag+1):n]
      cum_sums <- cumsum(sum_terms)
      gamma_ls <- cum_sums
      
      sigmas <- sigmas + 2*wkern(lag,lags)*gamma_ls[(lags+1-lag):length(gamma_ls)]
    }
    sigmas <- 1/hv * sigmas/kernels[(1+lags):n]^2 #rescaling
    
    return(sigmas)
  }
  
  sigmas <- foreach(path=1:paths, .combine = 'rbind')  %do% {sigma_func(path)}
  
  #Pad with zero in beginning for correct size
  if (paths == 1) {
    sigmas <- c(rep(0,lags),sigmas)
  } else {
    sigmas <- cbind(matrix(0,nrow = paths,ncol = lags),sigmas) 
  }
  
  #print(Sys.time()-p0)
  
  return(list(time = data$time[-1], sig = sigmas)) #Don't include time zero as this is not in dy
}

#Without parallel - faster that with parallel
est.mu.mat.2.0 <- function(data, hd, kern = kern.leftexp, wkern=kern.parzen){
  
  #hd <- 5 / (60*24*7*52)
  
  p0 <- Sys.time()
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hd)
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  mus <- matrix(nrow = paths, ncol = n)
  
  dy <- data$Y

  for (path in 1:paths) {
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy
    
    #zero lag
    sum_terms <- products
    mu_non_scaled <- 1/hd * cumsum(sum_terms)
    mus[path,] <- mu_non_scaled/kernels
  }

  print(Sys.time()-p0)
  
  return(list(time = data$time[-1],mu = mus)) #Don't include time zero, because dy doesn't include
}


