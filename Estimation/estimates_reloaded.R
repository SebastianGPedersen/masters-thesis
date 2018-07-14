#Install required packages if not already installed
packages <- c('foreach', 'doParallel')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load packages
invisible(sapply(packages, require, character.only = TRUE))

#Load sources
setwd(Sys.getenv("masters-thesis"))
source("kernels/kernels.R")
source("simulation/parameters.R") #gets number of logic units
registerDoParallel(n_logic_units)


######### FUNCTIONS ######### 

### SIGMA ###
#With parallel
est.sigma.mat.2.0 <- function(data, hv, kern = kern.leftexp, wkern=kern.parzen, lag=10, bandwidth_rescale = F){
  Heston <- data #possible fix?
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
  #rescaling <- kern((data$time[2:(length(data$time))]-t_now)/hv)
  rescaling <- kernels
  
  #Extra rescaling
  if (bandwidth_rescale) {
    K2 <- 0.5
    x <- (data$time[1]-data$time[2:length(data$time)])/hv
    scaling_integral <- 0.5 * (1-exp(2*x))
    
    rescaling <- rescaling / sqrt(0.5/scaling_integral) #sqrt because rescaling is squared later
  }
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  lags <- lag
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  #temp_func
  sigma_func <- function(path) {
    
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy
    
    # if (path == 1) {
    # print(dy[1:(lags+1)])
    # print(kernels[1:(lags+1)])
    # print(dy[1:(lags+1)]*kernels[1:(lags+1)])
    # print(sum((dy[1:(lags+1)]*kernels[1:(lags+1)])^2))
    # }
    
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
    sigmas <- 1/hv * sigmas/rescaling[(1+lags):n]^2 #rescaling
    
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
est.mu.mat.2.0 <- function(data, hd, kern = kern.leftexp, bandwidth_rescale = F){
  #data <- BS
  #hd <- hd_list[hd]
  
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")

  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hd)
  #rescaling <- kern((data$time[2:length(data$time)]-t_now)/hd)
  rescaling <- kernels
  
  #Extra rescaling
  if (bandwidth_rescale) {
    K2 <- 0.5
    x <- (data$time[1]-data$time[2:length(data$time)])/hd
    scaling_integral <- 0.5 * (1-exp(2*x))
    
    rescaling <- rescaling / sqrt(0.5/scaling_integral)
  }

  #Initialize for loop
  paths <- dim(data$Y)[1]
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  mus <- matrix(nrow = paths, ncol = n)

  for (path in 1:paths) {
    #path <- 1
    dy <- data$Y[path,]
    sum_terms <- kernels*dy
    
    mu_non_scaled <- 1/hd * cumsum(sum_terms)
    mus[path,] <- mu_non_scaled/rescaling
  }

  return(list(time = data$time[-1],mu = mus)) #Don't include time zero, because dy doesn't include
}

#Added a raw estimator as well
est.sigma.raw.mat.2.0 <- function(data, hv, kern = kern.leftexp){
  #data <- Heston
  #hd <- h_mu
  
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hv)
  #rescaling <- kern((data$time[2:length(data$time)]-t_now)/hd)
  rescaling <- kernels
  
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  sigmas <- matrix(nrow = paths, ncol = n)
  
  dy <- data$Y
  
  for (path in 1:paths) {
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy^2
    
    #zero lag
    sum_terms <- products
    sigmas_non_scaled <- 1/hv * cumsum(sum_terms)
    sigmas[path,] <- sigmas_non_scaled/rescaling
  }
  
  return(list(time = data$time[-1],sig = sigmas)) #Don't include time zero, because dy doesn't include
}

#Without parallel - faster that with parallel
est.mu.pre_avg.mat.2.0 <- function(data, hd, k_n, kern = kern.leftexp, bandwidth_rescale = F){
  #data <- path
  
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hd)
  #rescaling <- kern((data$time[2:length(data$time)]-t_now)/hd)
  rescaling <- kernels
  
  #Extra rescaling
  if (bandwidth_rescale) {
    K2 <- 0.5
    x <- (data$time[1]-data$time[2:length(data$time)])/hd
    scaling_integral <- 0.5 * (1-exp(2*x))
    
    rescaling <- rescaling / sqrt(0.5/scaling_integral)
  }
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  mus <- matrix(nrow = paths, ncol = n)
  
  for (path in 1:paths) {
    #path <- 1
    dy <- data$Y[path,]
    products <- kernels*dy
    pre_avg <- c(rep(0,k_n),est.NewPreAverage(delta_y = products, k_n)) #Pad with zeros in beggining so vector is corect size.
    sum_terms <- pre_avg
    
    mu_non_scaled <- 1/hd * cumsum(sum_terms)
    mus[path,] <- mu_non_scaled/rescaling
  }
  
  return(list(time = data$time[-1],mu = mus)) #Don't include time zero, because dy doesn't include
}

#Added a raw estimator as well
est.sigma.pre_avg.mat.2.0 <- function(data, hv, k_n, kern = kern.leftexp){
  #data <- Heston
  #hd <- h_mu
  
  #kern handling
  if(is.list(kern)) kern<-kern$kern
  if(!is.function(kern)) stop("kern should be either function or list containing function")
  
  t_now <- data$time[length(data$time)]
  
  #kernels
  kernels <- kern((data$time[1:(length(data$time)-1)]-t_now)/hv)
  #rescaling <- kern((data$time[2:length(data$time)]-t_now)/hd)
  rescaling <- kernels
  
  
  #Initialize for loop
  paths <- dim(data$Y)[1]
  
  n <- length(data$time)-1 #time includes time 0 and time t, wheras dy has one less
  
  sigmas <- matrix(nrow = paths, ncol = n)
  
  dy <- data$Y
  
  for (path in 1:paths) {
    #path <- 1
    dy <- data$Y[path,]
    products <- (kernels*dy)^2
    pre_avg <- c(rep(0,k_n),est.NewPreAverage(delta_y = products, k_n)) #Pad with zeros in beggining so vector is corect size.
    
    #zero lag
    sum_terms <- pre_avg
    sigmas_non_scaled <- 1/hv * cumsum(sum_terms)
    sigmas[path,] <- sigmas_non_scaled/rescaling^2
  }
  
  return(list(time = data$time[-1],sig = sigmas)) #Don't include time zero, because dy doesn't include
}


#This is just pseudo-code for the report
mu_estimator <- function(dy_vector, K_function, h_mu,time_points){
  
  t_end <- time_points[length(time_points)]
  kernels <- K_function(
              (time_points[1:(length(time_points)-1)]-t_end)/h_mu)
  sum_terms <- kernels*dy_vector
  mu_non_scaled <- 1/h_mu * cumsum(sum_terms)
  mu_estimates <- mu_non_scaled/kernels
  
  return(mu_estimates)
}
