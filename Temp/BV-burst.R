setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/BlackScholes.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("kernels/kernels.R")
source("vol/vol_estimators.R")

############################
# FACT OR FRICTION TABLE 2 #
############################

vol.est.performance <- function(data, thetas){
  # CALCULATES RV*/IV AND BV*/IV - FACT OR FRICTION TABLE 2
  estimation <- function(data, theta){
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
    return(data.table(theta = theta, IV = mean(IV), BV = mean(BVstar), BVIV = mean(BVstar/IV)))
  }
  # DO SO FOR EACH THETA
  mid <- NULL
  for(theta in thetas){
    mid <- rbind(mid, estimation(data, theta)) # naive append - shouldn't be bottleneck
  }
  return(mid)
}

thetas <- c(0.5, 1, 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HESTON-SIMULATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nsteps <- 100000; Npath <- 100 # 40 000; 10 000 # TAKES ROUGHLY 1.5 hours for all schemes

# MEMORY LOOP
total <- NULL
for(mem in 1:100){
  gc()
  print(mem)
  print(paste0("Loop cycle: ",mem," - ",Sys.time()))
  # BURST/JUMP SETTINGS
  burst_time <- 0.5; interval_length <- 0.05; alpha <- 0.8; beta <- 0.45;
  c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha); c_2 <- sqrt((1-2*beta)*0.001^2/(10/(60*24*7*52))^(1-2*beta));

  
  # RAW
  hest <- sim.heston(sim.setup(Nsteps = Nsteps, Npath = Npath))
  if(T){
    #JUMP
    hestJ<-sim.addjump(Heston_res = hest, burst_time = burst_time, interval_length = interval_length,
                       c_1 = c_1, alpha = alpha)
    # VB
    hestvb <- sim.addvb(Heston_res = hest, burst_time = burst_time, interval_length = interval_length,
                        c_2 = c_2, beta = beta, reverse = F, recenter = F)

    # VB DriftBurst - high alpha
    hestdbvb<-sim.adddb(Heston_res = hestvb, burst_time = burst_time, interval_length = interval_length,
                      c_1 = c_1, alpha = alpha, reverse = F)
    
    # VB DriftBurst - low alpha
    hestdbvb2<-sim.adddb(Heston_res = hestvb, burst_time = burst_time, interval_length = interval_length,
                        c_1 = (1-0.55)*0.005/(10/(60*24*7*52))^(1-0.55), alpha = 0.55, reverse = F)
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # TIME 2 GO
  
  # CALC AND TAG
  raw <- vol.est.performance(hest, thetas)
  raw <- raw[, id:= "raw"]
  
  J <- vol.est.performance(hestJ, thetas)
  J<-J[, id:= "J"]
  
  vb <- vol.est.performance(hestvb, thetas)
  vb<-vb[, id:= "vb"]
  
  dbvb <- vol.est.performance(hestdbvb, thetas)
  dbvb<-dbvb[, id:= "dbvb"]
  
  dbvb2 <- vol.est.performance(hestdbvb2, thetas)
  dbvb2<-dbvb2[, id:= "dbvb2"]

  total <- rbind(total, raw, J, vb, vbJ, dbvb, dbvb2)
}

# MEAN OF EACH
fullmean <- total[, lapply(.SD, mean), by = c("theta", "id"), .SDcols = c("BV")]

saveRDS(fullmean, file = "BV-burstOHNE-RV.rds")
