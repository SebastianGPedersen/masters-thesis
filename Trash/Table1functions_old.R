setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/estimates_revolution.R")
source("estimation/teststat.R")
source("estimation/rho.R")
source("estimation/rescaling.R")
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

# CREATE LIST OF FUN AND A LIST OF BURSTSETS

# This returns simulated raw | vb | db | vbdb | J |  vbJ  - all with the same burst settings
#funs <- list(raw, sim.addvb_arglist, sim.adddb_arglist, vbdb, sim.addjump_arglist, vbJ)
#args <- list(burstset, burstset, burstset, burstset, burstset, burstset)

# (ideally we would only give burst settings, and it would automatically provide the functions)
# (that way it is more intuitively similar to table1 - but it may not be that smart...)
# (what happens if we give burstsetting with beta = 0 and the VB uses that - might encounter weird shit like 0/0 etc.)
# (HERE we would have to control it outselves)

Table1 <- function(hest_setup, fun, args, h_list, ratio, t.index, lag, conf = 0.95){
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
  
  #sim <- list()
  # CREATE LIST OF SIMS
  #for(i in 1:length(funs)){
  #  sim <- append(sim, combine(hest, funs[[i]], args[[i]])) # this fucks up
  #}
  
  M <- length(args)
  out <- matrix(NA, nrow = M, ncol = length(h_list))
  colnames(out) <- paste0(hset*52*7*24*60*60)
  for(i in 1:M){
    print(i)
    print(Sys.time())
    sim <- combine(hest, funs[[i]], args[[i]])
    out[i,] <- Table1_estimation(sim_list = sim, h_list = h_list,
                                 ratio = ratio, t.index = t.index, lag = lag, conf = 0.95)
  }
  return(out)
}

# T ESTIMATION PART
Table1_estimation <- function(sim_list, h_list, ratio, t.index, lag, conf = 0.95){
  # HERE WE GO
  path <- sim_list
  # HIGHWAY
  # CLEAN SLATE
  N<-dim(path$Y)[1]
  Tstar <- numeric(N);  rho<-numeric(N); m <- numeric(N)
  
  #Transform Y to dY in path$Y
  path$Y <- t(diff(t(as.matrix(path$Y))))
  
  output <- numeric(length(h_list))
  ######## CALCULATE T estimator ##########
  for (h_index in 1:length(h_list)) {
    #p0 <- Sys.time()
    #print(paste0("memory = ",memory, ", path = ",j, ", ratio_index = ", ratio_index, sep = ""))
    
    mu_hat <- sqrt(h_list[h_index])*est.mu.mat.2.0(data = path, hd = h_list[h_index])$mu[,t.index]#,t.index = t.index)$mu
    sigma_hat2 <- est.sigma.mat.3.0(data = path, hv = ratio*h_list[h_index])$sig[,t.index]#,t.index = t.index,lag = lag)$sig
    
    # RESCALE
    mu <- est.rescale.mu(mu_vector = mu_hat, time_points =path$time[t.index], t_beginning = path$time[1], bandwidth = h_list[h_index])
    sigma2 <- est.rescale.sigma(sigma2_vector = sigma_hat2, time_points =path$time[t.index], t_beginning = path$time[1], bandwidth = ratio*h_list[h_index])
    
    Tstat <- mu/sqrt(sigma2)
    
    # Calculate T* for each subpath
    for(subpath in 1:N){
      Tstar[subpath] <- max(abs(Tstat[subpath,]))
      # fit rho
      rho[subpath] <- est.rho.MLE(Tstat[subpath,])
    }
    
    m <- rep(dim(Tstat)[2],N)
    
    # QUANTILES for all paths
    z<-est.z_quantile(m, rho, conf)$qZm
  
    #output <- list(output, mean(Tstar>=z))
    
    output[h_index] <- mean(Tstar>=z)
    #output <- mean(Tstar)
      
    #print((Sys.time()-p0)*length(ratio_list)*length(all_paths)*n_loops)
  }
  return(output)
}

