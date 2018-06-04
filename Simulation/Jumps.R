sim.addjump <- function(Heston_res, size, burst_time = 0.5, interval_length = 0.05, c_1 = 3, alpha = 0.75) {
  # Adds a jump corresponding to the given burst settings
  # OR simply jump of size: "size" at burst_time"
  
  #Intervals
  burst_begin_perc = burst_time-interval_length/2
  burst_end_perc = burst_time+interval_length/2
  
  burst_begin = burst_begin_perc * Heston_res$time[length(Heston_res$time)]
  burst_end = burst_end_perc * Heston_res$time[length(Heston_res$time)]
  tau = (burst_end+burst_begin)/2
  
  
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
  dt = Heston_res$time[2] - Heston_res$time[1]
  
  mu_add = vector(length=steps+1)
  mu_add[1] = 0
  
  #Calculate the mean vector (that has to be added to X and Y)
  for (i in 2:(steps+1)) {
    mu_add[i] = mu_add[i-1]+mu(Heston_res$time[i-1])*dt
  }
  # Jumpify
  ind<-which.max(abs(mu_add))
  
  if(missing(size)){
    size <- mu_add[ind]
  }
  mu_add[] <- 0
  mu_add[ind:length(mu_add)] <- size
  
  if ("X" %in% names(Heston_res)) {
    Heston_res$X = Heston_res$X+t(replicate(nrow(Heston_res$X), mu_add))
  }
  Heston_res$Y = Heston_res$Y+t(replicate(nrow(Heston_res$Y), mu_add))
  
  return(Heston_res)
}

# LIST ARG OVERLOAD
sim.addjump_arglist <- function(Heston_res, burstsetting) {
  # Adds a jump corresponding to the given burst settings
  # OR simply jump of size: "size" at burst_time"
  alpha <- burstsetting$alpha;  beta <- burstsetting$beta;
  burst_time <- burstsetting$burst_time;  interval_length <- burstsetting$interval_length;
  c_1 <- burstsetting$c_1;  c_2 <- burstsetting$c_2; reverse <- burstsetting$reverse; recenter <- burstsetting$recenter
  #Intervals
  burst_begin_perc = burst_time-interval_length/2
  burst_end_perc = burst_time+interval_length/2
  
  burst_begin = burst_begin_perc * Heston_res$time[length(Heston_res$time)]
  burst_end = burst_end_perc * Heston_res$time[length(Heston_res$time)]
  tau = (burst_end+burst_begin)/2
  
  
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
  dt = Heston_res$time[2] - Heston_res$time[1]
  
  mu_add = vector(length=steps+1)
  mu_add[1] = 0
  
  #Calculate the mean vector (that has to be added to X and Y)
  for (i in 2:(steps+1)) {
    mu_add[i] = mu_add[i-1]+mu(Heston_res$time[i-1])*dt
  }
  # Jumpify
  ind<-which.max(abs(mu_add))
  
  size <- mu_add[ind]
  
  mu_add[] <- 0
  mu_add[ind:length(mu_add)] <- size
  
  Heston_res$X = Heston_res$X+t(replicate(nrow(Heston_res$X), mu_add))
  Heston_res$Y = Heston_res$Y+t(replicate(nrow(Heston_res$X), mu_add))
  
  return(Heston_res)
}
