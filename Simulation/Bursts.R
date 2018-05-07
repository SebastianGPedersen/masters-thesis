source(paste(Sys.getenv("masters-thesis"),"Simulation/Heston.R",sep="/"))



# Nedenstående funktioner køres som følger:

# Heston_res = sim.heston(settings)
# Heston_vb = sim.addvb(Heston_res, burst_time, interval_length, c_2, beta)
# Heston_vbdb = sim.adddb(Heston_vb, burst_time, interval_length, c_1, alpha)

#OBS: burst_time og interval_length skal være i % and intervallængde.

sim.burstsetting <- function(alpha, beta, burst_time = 0.5, interval_length = 0.05,
                         c_1 = 3, c_2= 0.15){
  return(list(alpha = alpha, beta = beta, burst_time = burst_time,
              interval_length = interval_length, c_1 = c_1, c_2 = c_2))
}

sim.adddb <- function(Heston_res, burst_time = 0.5, interval_length = 0.5, c_1 = 3, alpha = 0.75, reverse = T) {
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
  
  Heston_res$X = Heston_res$X+t(replicate(nrow(Heston_res$X),mu_add))
  Heston_res$Y = Heston_res$Y+t(replicate(nrow(Heston_res$X),mu_add))
  
  #Check that we end up in same point
  if(reverse){
    if ((mu_add[length(mu_add)] > 10^(-6)) | (mu_add[length(mu_add)] < -10^(-6))) {
      stop("We don't end up in same point - something is wrong with the mu-vector")
    }
  }
  return(Heston_res)
}

sim.addvb <- function(Heston_res, burst_time = 0.5, interval_length = 0.5, c_2 = 0.15, beta = 0.4, reverse = T, recenter = T) {
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
  
  sigma_add = matrix(nrow = nrow(Heston_res$X) ,ncol=(steps+1))
  sigma_add[,1] = 0
  
  #This time we need the dW to calculate sigma^{vb}*dW
  #dW = (Heston_res$X[,2:(steps+1)]-Heston_res$X[,1:steps]) / (sqrt(Heston_res$vol[,1:steps])*sqrt(dt)) #Isolate dW in eq. from Heston-function
  dW = (Heston_res$X[,2:(steps+1)]-Heston_res$X[,1:steps]) / (sqrt(Heston_res$vol[,1:steps])) #Isolate dW in eq. from Heston-function
  
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
  
  vol_burst_vector = sapply(Heston_res$time, function(x) sigma(x)^2)
  new_vol = Heston_res$vol + t(replicate(nrow(Heston_res$X),vol_burst_vector))
  
  #new_epsilon = epsilon_u_vol*sqrt(new_vol) #Multiplicates matrixes entrance-wise
  
  #Change X,Y and vol
  Heston_res$X = Heston_res$X+sigma_add
  #Heston_res$Y = Heston_res$X+new_epsilon
  Heston_res$Y = Heston_res$Y+sigma_add
  Heston_res$vol = new_vol
  
  #Return
  return(Heston_res)
}