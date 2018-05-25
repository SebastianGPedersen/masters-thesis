setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("kernels/kernels.R")

# HELPER FUNCTIONS
raw <- function(hest, burstsetting){
  return(hest)
}

vbdb <- function(hest, burstsetting){
  mid <- sim.addvb_arglist(hest, burstsetting)
  final <- sim.adddb_arglist(mid, burstsetting)
  return(final)
}

vbJ <- function(hest, burstsetting){
  mid <- sim.addvb_arglist(hest, burstsetting)
  final <- sim.addjump_arglist(mid, burstsetting)
  return(final)
}

#Constants for bursts
c_1_func <- function(alpha) {(1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)}
c_2_func <- function(beta) {sqrt((1-2*beta)*0.001^2/(10/(60*24*7*52))^(1-2*beta))}

# </HELPER FUNCTIONS>

# SIMULATION PREP

#Burst settings1
alpha <- 0.8
beta <- 0.1
burst_time = 0.5; interval_length = 0.05



# CREATE BURST SETTING
burstset <- sim.burstsetting(alpha = alpha, beta = beta,
                             burst_time = burst_time, interval_length = interval_length,
                             c_1 = c_1_func(alpha), c_2 = c_2_func(beta))

# CREATE LIST OF FUN AND A LIST OF BURSTSETS

# This returns simulated raw | vb | db | vbdb | J |  vbJ  - all with the same burst settings
funs <- list(raw, sim.addvb_arglist, sim.adddb_arglist, vbdb, sim.addjump_arglist, vbJ)
args <- list(burstset, burstset, burstset, burstset, burstset, burstset)

# (ideally we would only give burst settings, and it would automatically provide the functions)
# (that way it is more intuitively similar to table1 - but it may not be that smart...)
# (what happens if we give burstsetting with beta = 0 and the VB uses that - might encounter weird shit like 0/0 etc.)
# (HERE we would have to control it outselves)

Table1_simulation <- function(hest_setup, fun, args){
  # WE NEED TO START WITH A RAW
  hest <- sim.heston(hest_setup)
  
  # COMBINE FUNCTION AND SETTINGS
  combine <- function(hest, fun, burstsetting){
    # DATA SHOULD BE RAW HESTON
    # FUN COULD BE Add.vbdb
    # PARAMS IS A BURSTSET LIST
    mid <- fun(hest, burstsetting)
    # Remove X/Vol
    mid$X <- NULL
    mid$vol <- NULL
    return(mid)
  }
  
  sim <- list()
  # CREATE LIST OF SIMS
  for(i in 1:length(funs)){
    sim <- list(sim, combine(hest, funs[[i]], args[[i]]))
  }
  return(sim)
}


# T ESTIMATION PART
Table1_estimation <- function(sim_list, h_list, ratio, t.index, lag, conf = 0.95){
  # HERE WE GO
  all_paths <- sim_list
  # HIGHWAY
  
  # output list
  output <- list()
  for (j in 1:length(all_paths)) {
    #j <- 1
    path <- all_paths[[j]]
    
    # CLEAN SLATE
    N<-dim(path$time)[1]
    Tstar <- numeric(N);  rhom <- numeric(N);  rhorho<-numeric(N);
    
    #Transform Y to dY in path$Y
    path$Y <- t(diff(t(as.matrix(path$Y))))
    
    ######## CALCULATE T estimator ##########
    for (h_index in 1:length(h_list)) {
      p0 <- Sys.time()
      #print(paste0("memory = ",memory, ", path = ",j, ", ratio_index = ", ratio_index, sep = ""))
      
      mu_hat <- est.mu.mat(data = path, hd = h_list[h_index],t.index = t.index)$mu
      sigma_hat2 <- est.sigma.mat(data = path, hv = ratio*h_list[h_index],t.index = t.index,lag = lag)$sig
      Tstat <- abs(sqrt(h_list[h_index])*mu_hat/sqrt(sigma_hat2))
      
      # Calculate T* for each subpath
      for(subpath in 1:N){
        Tstar[subpath] <- max(abs(Tstat[subpath,]))
        # fit rho
        rho <- est.rho(Tstat$test[subpath,])
        rhom[subpath] <- rho$m
        rhorho[subpath] <- rho$rho
      }
      
      # QUANTILES for all paths
      z<-est.z_quantile(rhom, rhorho, conf)$qZm
      
      output <- list(output, mean(Tstar>=z))
      
      print(mean(Tstar>=z))
      
      #print((Sys.time()-p0)*length(ratio_list)*length(all_paths)*n_loops)
    }
    return(output)
  }
}
