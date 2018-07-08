#Install required packages if not already installed
packages <- c('foreach', 'Rcpp', 'doParallel')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load packages
invisible(sapply(packages, require, character.only = TRUE))

#Load sources
setwd(Sys.getenv("masters-thesis"))
source("kernels/kernels.R")
sourceCpp('Estimation/sigma.cpp')
source("simulation/parameters.R") #gets number of logic units
registerDoParallel(n_logic_units)

######### FUNCTIONS ######### 

### SIGMA ###
#With parallel
est.sigma.mat.3.0 <- function(data, hv, kern = kern.leftexp, wkern=kern.parzen, lag=10, bandwidth_rescale = F){
  #data <- path
  #hv <- 5 / (60*24*7*52)
  #lag <- 10
  
  #p0 <- Sys.time()
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  if(is.list(wkern)) wkern<-wkern$kern
  if(!is.function(wkern)) stop("wkern should be either function or list containing function")
  
  
  #For sigma3.0
  t_now <- data$time[length(data$time)]
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hv)
  #rescaling <- kern((data$time[2:length(data$time)]-t_now)/hv)
  rescaling <- kernels
  
  #Extra rescaling
  if (bandwidth_rescale) {
    K2 <- 0.5
    x <- (data$time[1]-data$time[2:length(data$time)])/hv
    scaling_integral <- 0.5 * (1-exp(2*x))
    
    rescaling <- rescaling / sqrt(0.5/scaling_integral) #sqrt because it is squared later
  }
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  lags <- lag
  
  #temp_func
  sigma_func <- function(path) {
    
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy
    sigmas_non_scaled <- sigmas_cpp(KdY = products,lags = lags)
    #return(sigmas_non_scaled)
    sigmas <- 1/hv * sigmas_non_scaled/rescaling^2 #rescaling
    return(sigmas)
  }
  
  sigmas <- foreach(path=1:paths, .combine = 'rbind')  %dopar% {sigma_func(path)}
  
  #for(path in 1:paths) {sigmas[path,] <- sigma_func(path)}
  
  return(list(time = data$time[-1], sig = sigmas)) #Don't include time zero as this is not in dy
}


#Pseudo-code for report
sigma_estimator <- function(dY_vector, h_sigma, K_func, 
                            time_points,lags=10){

  
  t_end <- time_points[length(time_points)]
  kernels <- K_func(
              (time_points[1:(length(time_points)-1)]-t_end)/h_sigma)
  products <- kernels*dy_vector
  
  #Using the 'sigmas_cpp' function
  sigmas_non_scaled <- sigmas_cpp(KdY = products,lags = lags)
  sigmas <- 1/hv * sigmas_non_scaled/kernels^2
  return(sigmas)
}



