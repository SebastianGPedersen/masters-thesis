#Install required packages if not already installed
packages <- c('foreach', 'Rcpp')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load packages
invisible(sapply(packages, require, character.only = TRUE))

#Load sources
setwd(Sys.getenv("masters-thesis"))
source("kernels/kernels.R")
sourceCpp('Estimation/sigma.cpp')


######### FUNCTIONS ######### 

### SIGMA ###
#With parallel
est.sigma.mat.3.0 <- function(data, hv, kern = kern.leftexp, wkern=kern.parzen, lag=100){
  #data <- path
  #hv <- 5 / (60*24*7*52)
  #lag <- 100
  
  #p0 <- Sys.time()
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  
  #For sigma3.0
  kernels <- kern((data$time[1:(length(data$time)-1)]-data$time[2:length(data$time)])/hv)
  
  #Initialize for loop
  paths <- dim(data$Y)[1]

  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  #sigmas <- matrix(nrow = paths, ncol = n)
  
  #temp_func
  sigma_func <- function(path) {
    
    #path <- 1
    dy <- data$Y[path,]
    return(sigmas_cpp(dY = dy,kernels = kernels,lags = lag, bandwidth = hv))
  }
  
  sigmas <- foreach(path=1:paths, .combine = 'rbind')  %do% {sigma_func(path)}
  
  #for(path in 1:paths) {sigmas[path,] <- sigma_func(path)}
  
  return(list(time = data$time[-1], sig = sigmas)) #Don't include time zero as this is not in dy
}




