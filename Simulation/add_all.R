#Load sources
setwd(Sys.getenv("masters-thesis"))
source("simulation/bursts.R")
source("simulation/jumps.R")
library(mgcv) #Used to extract unique rows from matrix

#Add the volatility bursts first
sim.add_all <- function(Heston, burst_args) {
  #burst_args <- burstsettings
  
  
  #Initialize result_list
  all_sims <- list()
  for (i in 1:length(burst_args)) {
    all_sims[[i]] <- list(Y = Heston$Y, X = Heston$X, vol = Heston$vol, time = Heston$time, params = burst_args[i])
  }

  ### Add volatility - this MUST be done before drift burst and jump (this has SIGNIFICANT computation time)
  #Get all beta's
  betas_c2s <- matrix(0,nrow = length(burst_args),ncol = 4)
  for (i in 1:length(burst_args)) {
    betas_c2s[i,1] <- burst_args[[i]]$beta
    betas_c2s[i,2] <- burst_args[[i]]$c_2
    betas_c2s[i,3] <- burst_args[[i]]$recenter
    betas_c2s[i,4] <- burst_args[[i]]$reverse
  }
  
  #Get unique combinations of c(beta,c2)
  betas_c2s_unique <- uniquecombs(betas_c2s)
  if (!is.null(nrow(betas_c2s_unique))) {
    betas_c2s_u_zero <- betas_c2s_unique[betas_c2s_unique[,1] != 0,] #Remove the zero.
    
    if (is.null(nrow(betas_c2s_u_zero))) {
      betas_c2s_u_zero <- matrix(betas_c2s_u_zero, ncol = 4) #If only one row -> create matrix
    }
    
    #Add all combinations of c(beta,c_2) to Heston
    for (row in 1:nrow(betas_c2s_u_zero)) {
      Heston_vb <- sim.addvb.2.0(Heston_res = Heston, 
                                 beta = betas_c2s_u_zero[row,1], 
                                 c_2 = betas_c2s_u_zero[row,2],
                                 recenter = betas_c2s_u_zero[row,3],
                                 reverse = betas_c2s_u_zero[row,4])
      for (replace in 1:nrow(betas_c2s)) {
        if ((betas_c2s[replace,1] ==  betas_c2s_u_zero[row,1]) & (betas_c2s[replace,2] ==  betas_c2s_u_zero[row,2])) {
          all_sims[[replace]]$Y <- Heston_vb$Y
        }
      }
    }
  }
  
  ### Add drift bursts OR jump to all (this has not insignificant computation time, so should be optimized if correct)
  for (i in 1:length(burst_args)) {
    #i <- 1
    alpha <- burst_args[[i]]$alpha
    jump <- burst_args[[i]]$jump
    c_1 <- burst_args[[i]]$c_1
    
    if (alpha > 0 & jump) { #add 
      all_sims[[i]]$Y <- sim.addjump(Heston_res = all_sims[[i]],alpha = alpha, c_1 = c_1)$Y
    } else if (alpha > 0) { #add drift burst
      all_sims[[i]]$Y <- sim.adddb(Heston_res = all_sims[[i]],alpha = alpha, c_1 = c_1, reverse = burst_args[[i]]$reverse)$Y    
    }
  }

  #return
  return(all_sims)
}

#Add everything to every burst
sim.add_all2 <- function(Heston, burst_args) {
  #burst_args <- burstsettings
  
  
  #Initialize result_list
  all_sims <- list()
  for (i in 1:length(burst_args)) {
    all_sims[[i]] <- list(Y = Heston$Y, X = Heston$X, vol = Heston$vol, time = Heston$time, params = burst_args[i])
  }

  
  ### Add the bursts
  for (i in 1:length(burst_args)) {
    #i <- 
    beta <- burst_args[[i]]$beta
    alpha <- burst_args[[i]]$alpha
    jump <- burst_args[[i]]$jump
    
    if (beta > 0) {
      all_sims[[i]]$Y <- sim.addvb.2.0(Heston_res = all_sims[[i]],
                                       beta = beta, 
                                       c_2 = burst_args[[i]]$c_2,
                                       burst_time = burst_args[[i]]$burst_time,
                                       reverse = burst_args[[i]]$reverse,
                                       recenter = burst_args[[i]]$recenter)$Y
    }
    
    if (alpha > 0 & jump) { #add 
      all_sims[[i]]$Y <- sim.addjump(Heston_res = all_sims[[i]],
                                     alpha = alpha, 
                                     c_1 = burst_args[[i]]$c_1,
                                     burst_time = burst_args[[i]]$burst_time)$Y
      
    } else if (alpha > 0) { #add drift burst
      all_sims[[i]]$Y <- sim.adddb(Heston_res = all_sims[[i]],
                                   alpha = alpha, 
                                   c_1 = burst_args[[i]]$c_1,
                                   burst_time = burst_args[[i]]$burst_time,
                                   reverse = burst_args[[i]]$reverse)$Y    
    }
  }
  
  #return
  return(all_sims)
}