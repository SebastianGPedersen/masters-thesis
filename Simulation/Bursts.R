#Install required packages if not already installed
packages <- c('foreach', 'doParallel', 'Rcpp')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load packages
invisible(sapply(packages, require, character.only = TRUE))

#Load sources
setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/parameters.R") #gets number of logic units
sourceCpp("simulation/addvb.cpp")
registerDoParallel(n_logic_units)


# Nedenstående funktioner køres som følger:

# Heston_res = sim.heston(settings)
# Heston_vb = sim.addvb(Heston_res, burst_time, interval_length, c_2, beta)
# Heston_vbdb = sim.adddb(Heston_vb, burst_time, interval_length, c_1, alpha)

#OBS: burst_time og interval_length skal være i % and intervallængde.

sim.burstsetting <- function(alpha, beta, jump = F, burst_time = 0.5, interval_length = 0.05,
                         c_1 = 3, c_2= 0.15, reverse = F, recenter = F){
  return(list(jump = jump, alpha = alpha, beta = beta, burst_time = burst_time,
              interval_length = interval_length, c_1 = c_1, c_2 = c_2, reverse = reverse, recenter = recenter))
}

sim.adddb <- function(Heston_res, burst_time = 0.5, interval_length = 0.05, c_1 = 3, alpha = 0.75, reverse = F) {
  #Intervals
  burst_begin_perc = burst_time-interval_length/2
  if(reverse){
    burst_end_perc = burst_time+interval_length/2
  }
  else{
    burst_end_perc = burst_time
  }
  
  burst_begin = burst_begin_perc * Heston_res$time[length(Heston_res$time)]
  burst_end = burst_end_perc * Heston_res$time[length(Heston_res$time)]
  
  if(reverse){
    tau = (burst_end+burst_begin)/2
  }
  else{
    tau = burst_end
  }
  
  #Define mu-function
  mu <- function(t) {
    if ((t>=burst_begin) & (t <=burst_end) & (t !=tau)){
      mu_t = c_1*sign(t-tau)/abs(tau-t)^alpha
    } else {
      mu_t=0
    }
    return(mu_t)
  }
  
  #initializations
  steps = length(Heston_res$time)-1
  #dt = Heston_res$time[2] - Heston_res$time[1] # this should be generalized if we do uneven dt
  dt = diff(Heston_res$time)
  
  mu_add = vector(length=steps+1)
  mu_add[1] = 0
  
  #Calculate the mean vector (that has to be added to X and Y)
  for (i in 2:(steps+1)) {
    mu_add[i] = mu_add[i-1]+mu(Heston_res$time[i-1])*dt[i-1]
  }
  
  if ("X" %in% names(Heston_res)){
    Heston_res$X <- Heston_res$X+t(replicate(nrow(Heston_res$X),mu_add))
    
  }
  Heston_res$Y = Heston_res$Y+t(replicate(nrow(Heston_res$Y),mu_add))
  
  #Check that we end up in same point
  if(reverse){
    if ((mu_add[length(mu_add)] > 10^(-6)) | (mu_add[length(mu_add)] < -10^(-6))) {
      stop("We don't end up in same point - something is wrong with the mu-vector")
    }
  }
  return(Heston_res)
}

sim.addvb <- function(Heston_res, burst_time = 0.5, interval_length = 0.05, c_2 = 0.15, beta = 0.4, reverse = F, recenter = F) {
  #Intervals
  burst_begin_perc = burst_time-interval_length/2
  if(reverse){
    burst_end_perc = burst_time+interval_length/2
  }
  else{
    burst_end_perc = burst_time
  }
  burst_begin = burst_begin_perc * Heston_res$time[length(Heston_res$time)]
  burst_end = burst_end_perc * Heston_res$time[length(Heston_res$time)]
  if(reverse){
    tau = (burst_end+burst_begin)/2
  }
  else{
    tau = burst_end
  }
  
  #Define sigma-function
  sigma <- function(t) {
    if ((t>=burst_begin) & (t <=burst_end) & (t !=tau)){
        sigma_t = c_2/abs(tau-t)^beta
      } else {
        sigma_t = 0
    }
    return(sigma_t)
  }
  
  #initializations
  steps = length(Heston_res$time)-1
  #dt = Heston_res$time[2]-Heston_res$time[1]
  #dt = diff(Heston_res$time)
  paths <- nrow(Heston_res$X)
  sigma_add = matrix(nrow = paths ,ncol=(steps+1))
  sigma_add[,1] = 0
  
  #This time we need the dW to calculate sigma^{vb}*dW
  #dW = (Heston_res$X[,2:(steps+1)]-Heston_res$X[,1:steps]) / (sqrt(Heston_res$vol[,1:steps])*sqrt(dt)) #Isolate dW in eq. from Heston-function
  dW <- (Heston_res$X[,2:(steps+1)]-Heston_res$X[,1:steps]) / (sqrt(Heston_res$vol[,1:steps])) #Isolate dW in eq. from Heston-function

  #Calculate sum(sigma*dW) # KAN sqrt(dt) fjernes? De g�r ud med hinanden i sigma_add...
  for (i in 2:(steps+1)) {
    #sigma_add[,i] = sigma_add[,i-1]+sigma(Heston_res$time[i-1])*sqrt(dt[i-1])*dW[,i-1]
    sigma_add[,i] = sigma_add[,i-1]+sigma(Heston_res$time[i-1])*dW[,i-1]
  }

  #From footnote 11 we need to recenter the sigma_add so it reverts
  if (recenter) {
    index_vector = ((Heston_res$time >=burst_begin) & (Heston_res$time <=burst_end)) #get true-false vector of interval
    len = sum(index_vector)
    extract = sigma_add[,ncol(sigma_add)]/len          #this has to be extracted to recenter (different in every path, it is a vector)
    
    for (i in 1:(ncol(sigma_add)-1)){                  #if we are in interval, the sigma AFTER has to be altered
      sigma_add[,i+1] = sigma_add[,i+1] - sum(index_vector[1:i])*extract
    }
      
    #Now that sigma is recentered we check that every path is centered correct
    if ((max(sigma_add[,ncol(sigma_add)]) > 10^(-6)) | (max(sigma_add[,ncol(sigma_add)]) < -10^(-6))) {
      stop("We don't end up in same point - something is wrong with the sigma-vector")
    }
  }
    
  #We are changing the vol, so epsilon also changes
  #epsilon_u_vol = (Heston_res$Y-Heston_res$X)/sqrt(Heston_res$vol) #this is gamma/sqrt(n) * rnorm
  
  vol_burst_vector <- sapply(Heston_res$time, function(x) sigma(x)^2)
  new_vol <- Heston_res$vol + matrix(rep(vol_burst_vector,nrow(Heston_res$X)),nrow = paths, byrow = T) #FS 11-07 | byrow = T

  #new_epsilon = epsilon_u_vol*sqrt(new_vol) #Multiplicates matrixes entrance-wise
  
  #Change X,Y and vol
  Heston_res$X = Heston_res$X+sigma_add
  #Heston_res$Y = Heston_res$X+new_epsilon
  Heston_res$Y = Heston_res$Y+sigma_add
  Heston_res$vol = new_vol
  
  #Return
  return(Heston_res)
}
sim.addvb.2.0 <- function(Heston_res, burst_time = 0.5, interval_length = 0.05, c_2 = 0.15, beta = 0.4, reverse = F, recenter = F) {
  
  #Intervals
  burst_begin_perc = burst_time-interval_length/2
  if(reverse){
    burst_end_perc = burst_time+interval_length/2
  }
  else{
    burst_end_perc = burst_time
  }
  burst_begin = burst_begin_perc * Heston_res$time[length(Heston_res$time)]
  burst_end = burst_end_perc * Heston_res$time[length(Heston_res$time)]
  if(reverse){
    tau = (burst_end+burst_begin)/2
  }
  else{
    tau = burst_end
  }
  
  #initializations
  steps <- length(Heston_res$time)-1
  paths <- nrow(Heston_res$X)
  
  #sigma_add <- matrix(nrow = paths ,ncol=(steps+1))

  #Get the sigma_vb to add
   t <- Heston_res$time[1:steps]
   sigma_t <- numeric(length = steps)
   indexes <- ((t>=burst_begin) & (t <=burst_end) & (t!=tau))
   sigma_t[indexes] <- c_2/abs(tau-t[indexes])^beta
  
  #This time we need the dW to calculate sigma^{vb}*dW
  # for (i in 1:paths) {
  #   dW <- (Heston_res$X[i,2:(steps+1)]-Heston_res$X[i,1:steps]) / (sqrt(Heston_res$vol[i,1:steps])) #Isolate dW in eq. from Heston-function
  #   temp <- sigma_t*dW
  #   sigma_add[i,] <- c(0,cumsum(temp))
  # }
  
  p0 <- Sys.time()
  
  sigma_add <- vol_add_cpp(X = Heston_res$X, vol = Heston_res$vol, sigma_t = sigma_t)
  
  # dW <- (Heston_res$X[,2:(steps+1)]-Heston_res$X[,1:steps]) / (sqrt(Heston_res$vol[,1:steps])) #Isolate dW in eq. from Heston-function
  # temp <- t(t(dW)*sigma_t)
  # sigma_add <- cbind(rep(0,paths),t(apply(temp,1,cumsum)))
  
  # path_func <- function(i) {
  #   return(c(0,cumsum(sigma_t*dW[i,])))
  # }
  # 
  # sigma_add <- foreach(path=1:paths, .combine = 'rbind')  %dopar% {path_func(path)}

  #From footnote 11 we need to recenter the sigma_add so it reverts
  if (recenter) {
    index_vector = ((Heston_res$time >=burst_begin) & (Heston_res$time <=burst_end)) #get true-false vector of interval
    len = sum(index_vector)
    extract = sigma_add[,ncol(sigma_add)]/len          #this has to be extracted to recenter (different in every path, it is a vector)
    
    for (i in 1:(ncol(sigma_add)-1)){                  #if we are in interval, the sigma AFTER has to be altered
      sigma_add[,i+1] = sigma_add[,i+1] - sum(index_vector[1:i])*extract
    }
    
    #Now that sigma is recentered we check that every path is centered correct
    if ((max(sigma_add[,ncol(sigma_add)]) > 10^(-6)) | (max(sigma_add[,ncol(sigma_add)]) < -10^(-6))) {
      stop("We don't end up in same point - something is wrong with the sigma-vector")
    }
  }
  
  #We are changing the vol, so epsilon also changes
  #epsilon_u_vol = (Heston_res$Y-Heston_res$X)/sqrt(Heston_res$vol) #this is gamma/sqrt(n) * rnorm
  
  #Get the sigma_vb to add
  vol_burst_vector <- numeric(length(Heston_res$time))
  t <- Heston_res$time
  sigma_t <- numeric(length = steps)
  indexes <- ((t>=burst_begin) & (t <=burst_end) & (t!=tau))
  vol_burst_vector[indexes] <- (c_2/abs(tau-t[indexes])^beta)^2

  new_vol <- Heston_res$vol + matrix(rep(vol_burst_vector,nrow(Heston_res$X)),nrow = paths, byrow = T) #FS 11-07 | byrow = T
  
  #new_epsilon = epsilon_u_vol*sqrt(new_vol) #Multiplicates matrixes entrance-wise
  
  #Change X,Y and vol
  if ("X" %in% names(Heston_res)) {
    Heston_res$X <- Heston_res$X+sigma_add
  }
  #Heston_res$Y = Heston_res$X+new_epsilon
  Heston_res$Y <- Heston_res$Y+sigma_add
  Heston_res$vol <- new_vol
  
  #Return
  return(Heston_res)
}

# Pass arguments as list
sim.adddb_arglist <- function(Heston, burstsetting) {
  # UNPACK
  alpha <- burstsetting$alpha;  beta <- burstsetting$beta;
  burst_time <- burstsetting$burst_time;  interval_length <- burstsetting$interval_length;
  c_1 <- burstsetting$c_1;  c_2 <- burstsetting$c_2; reverse <- burstsetting$reverse; recenter <- burstsetting$recenter
  
  Heston_res <- sim.adddb(Heston, burst_time = burst_time, interval_length = interval_length,
                          c_1 = c_1, alpha = alpha, reverse = reverse)
  
  return(Heston_res)
}

sim.addvb_arglist <- function(Heston, burstsetting) {
  #UNPACK
  # UNPACK
  alpha <- burstsetting$alpha;  beta <- burstsetting$beta;
  burst_time <- burstsetting$burst_time;  interval_length <- burstsetting$interval_length;
  c_1 <- burstsetting$c_1;  c_2 <- burstsetting$c_2; reverse <- burstsetting$reverse; recenter <- burstsetting$recenter
  
  Heston_res <- sim.addvb.2.0(Heston, burst_time = burst_time, interval_length = interval_length,
                          c_2 = c_2, beta = beta, reverse = reverse, recenter = recenter)
  
  #Return
  return(Heston_res)
}