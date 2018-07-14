setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/BlackScholes.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("kernels/kernels.R")
source("vol/vol_estimators.R")

############################
#         BV TABLE 1       #
############################

vol.est.performance <- function(hestdata, alphas, betas, c1s, c2s, thetas){
  # CALCULATES BV*/IV - FACT OR FRICTION TABLE 2
  estimation <- function(data, theta, alpha, beta){
    # SETUP ESTIMATION/CORE CALCULATIONS
    K <- function(theta, N){
      return(floor(theta*sqrt(N))+floor(theta*sqrt(N))%%2)
    }
    # IV
    dt <- data$time[2]-data$time[1]
    IV <- rowSums(data$vol*dt)
    
    # REAL PROCESS Y #
    Y <- data$Y
    N <- length(Y[1,])
    k <- K(theta, N)
    # PRE AVERAGING
    Rstar    <- matrix(NA, nrow = dim(Y)[1], ncol = N-k+1)
    omegaest <- rep(NA, dim(Y)[1]) #matrix(NA, nrow = dim(Y)[1], ncol = dim(Y)[2])
    BVstar   <- rep(NA, dim(Y)[1]) #matrix(NA, nrow = dim(Y)[1], ncol = dim(Y)[2])
    
    for(i in 1:dim(Y)[1]){
      Rstar[i,]    <- vol_est_preA(Y[i,], k)
      # OMEGA
      omegaest[i] <- vol.est.omega2(Data_subset = diff(Y[i,]))
      # BV*
      BVstar[i]   <- vol_est_BVstar(Rstar[i,], k, omegaest[i], theta)
    }
    return(data.table(theta = theta, IV = mean(IV), BV = mean(BVstar), BVIV = mean(BVstar/IV), alpha = alpha, beta = beta))
  }
  # DO SO FOR EACH MODE
  res <- NULL
  
  for(i in 1:length(betas)){
    if(betas[i] == 0){
      betadata <- hestdata
    }
    else{
      betadata <- sim.addvb(Heston_res = hestdata, burst_time = burst_time, interval_length = interval_length,
                          c_2 = c2s[i], beta = betas[i], reverse = F, recenter = F)
    }
  
    for(j in 1:length(alphas)){
      if(alphas[j] == 0){
        alphadata <- betadata
      }
      else{
        alphadata <- sim.adddb(Heston_res = betadata, burst_time = burst_time, interval_length = interval_length,
                          c_1 = c1s[j], alpha = alphas[j], reverse = F)
      }
      res <- rbind(res, estimation(alphadata, theta, alphas[j], betas[i])) # naive append - shouldn't be bottleneck
      gc()
    }
  }
  

  return(res)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nsteps <- 100000; Npath <- 100 # 40 000; 10 000 # TAKES ROUGHLY 1.5 hours for all schemes
theta <- 1

# BURST/JUMP SETTINGS
burst_time <- 0.5; interval_length <- 0.05;
# BURST SETTINGS (a = 0, 0.55, 0.65, 0.75 | b = 0.0, 0.1, 0.2, 0.3, 0.4), c_1 and c_2 from Batman
c_1_func <- function(alpha) {
  c_1 <- 0 #If alpha == 0
  if (alpha == 0.55) {
    c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)
  } 
  else if (alpha == 0.65){
    c_1 <- (1-alpha)*0.01/(10/(60*24*7*52))^(1-alpha)
  } 
  else if (alpha == 0.75) {
    c_1 <- (1-alpha)*0.015/(10/(60*24*7*52))^(1-alpha)
  }
  return(c_1)
}

c_2_func <- function(beta) {
  c_2 <- 0 #if beta = 0
  if (beta == 0.1) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.25)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.2) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.5)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.3) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.75)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.4) {
    c_2 <- sqrt((1-2*beta)*(0.00093*1)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  return(c_2)
}

alphas <- c(0, 0.55, 0.65, 0.75)
betas <-  c(0,0.1,0.2,0.3,0.4)
c1s    <- c(c_1_func(0), c_1_func(0.55), c_1_func(0.65), c_1_func(0.75) )
c2s    <- c(c_2_func(0), c_2_func(0.1), c_2_func(0.2), c_2_func(0.3), c_2_func(0.4) )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LET'S GO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MEMORY LOOP
total <- NULL
for(mem in 1:100){
  gc()
  print(mem)
  print(paste0("Loop cycle: ",mem," - ",Sys.time()))
  
  # RAW
  hest <- sim.heston(sim.setup(Nsteps = Nsteps, Npath = Npath))
  raw <- vol.est.performance(hest, alphas, betas, c1s, c2s, theta)
  total <- rbind(total, raw)
}

# MEAN OF EACH
fullmean <- total[, lapply(.SD, mean), by = c("alpha", "beta"), .SDcols = c("BV", "IV", "BVIV")]

saveRDS(fullmean, file = "BV-TABLE1.rds")
