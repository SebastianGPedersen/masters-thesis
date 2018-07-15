setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/estimates_reloaded.R")
source("estimation/estimates_revolution.R")
source("estimation/teststat.R")
source("estimation/rho.R")
source("estimation/rescaling.R")
source("kernels/kernels.R")
library(mgcv) #Used to extract unique rows from matrix

# T ESTIMATION PART
Table1_estimation <- function(sim_list, h_list, ratio, t.index, lag, conf = 0.95){
  #sim_list <- all_sims[[1]]
  # HERE WE GO
  
  #restructure to vector if not
  if(length(ratio) != length(h_list)) {
    ratio <- rep(ratio[1],length(h_list))
  }
  
  path <- sim_list
  # HIGHWAY
  # CLEAN SLATE
  N<-dim(path$Y)[1]
  Tstar <- numeric(N);  rho<-numeric(N); m <- numeric(N)
  test <- 1
  #Transform Y to dY in path$Y
  path$Y <- t(diff(t(as.matrix(path$Y))))
  
  output <- numeric(length(h_list))
  ######## CALCULATE T estimator ##########
  for (h_index in 1:length(h_list)) {
    #h_index <- 1
    #p0 <- Sys.time()
    #print(paste0("memory = ",memory, ", path = ",j, ", ratio_index = ", ratio_index, sep = ""))
    
    #2.0 and 3.0 can now also rescale bandwidth (almost twice as fast as seperate functions)
    mu <- sqrt(h_list[h_index])*est.mu.mat.2.0(data = path, hd = h_list[h_index], bandwidth_rescale = T)$mu[,t.index]#,t.index = t.index)$mu
    sigma2 <- est.sigma.mat.2.0(data = path, hv = ratio[h_index]*h_list[h_index], lag = lag, bandwidth_rescale = T)$sig[,t.index]#,t.index = t.index,lag = lag)$sig
    
    Tstat <- mu/sqrt(sigma2)
    
    #Calculate T* for each subpath
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

#This creates the structure of the table
Table1_func <- function(all_sims, h_list, ratio, t.index, lag, conf = 0.95){
  
  #all_sims <- all_simulations
  #Initialization
  M <- length(all_sims)
  out <- as.data.frame(matrix(nrow = M,ncol = 3+length(h_list)))
  
  #Create colnames and argument names
  h_colnames <- numeric(length(h_list))
  for (h in 1:length(h_list)){
    h_colnames[h] <- paste("h_n =",round(h_list[h]*(52*7*24*60),0),"min")
  }
  colnames(out) <- c("alpha","beta","jump",h_colnames)
  
  for (arg in 1:length(all_sims)) {
    #arg <- 1
    out[arg,1] <- all_sims[[arg]]$params[[1]]$alpha
    out[arg,2] <- all_sims[[arg]]$params[[1]]$beta
    out[arg,3] <- all_sims[[arg]]$params[[1]]$jump
  }
  
  #Add data
  for (i in 1:M) {
    out[i,3+1:length(h_list)] <- Table1_estimation(sim_list = all_sims[[i]], h_list = h_list,
                                 ratio = ratio, t.index = t.index, lag = lag, conf = 0.95)*100
  }
  
  return(out)
}

#Table restructuring
restructure_table1 <-function(table1, h_list) {
  #table1 <- Table1_results

  #Calculate number of rows (in output table)
  betas <- sort(unique(table1[,'beta']))
  n_rows <- length(betas)*length(h_list)
  
  #Calculate number of columns (in output table)
  temp <- table1[table1[,"alpha"] > 0,c("alpha","jump")]
  n_bursts <- uniquecombs(temp)
  n_cols <- 3+nrow(n_bursts) #bandwidth,beta,heston + n_bursts
  
  
  
  #### Create skeleton of output table
  output_table <- as.data.frame(matrix(nrow = n_rows, ncol = n_cols))
  #Bandwidths
  h_list <- sort(h_list)
  h_mins <- numeric(length(h_list))
  for (h in 1:length(h_list)){
    h_mins[h] <- paste("h_n =",round(h_list[h]*(52*7*24*60),0),"min")
  }
  output_table[,1] <- rep(h_mins,each = length(betas))
  #Betas
  output_table[,2] <- rep(betas,length(h_list))
  #Columns names (burst processes)
  db <- n_bursts[(n_bursts[,1] > 0) & (n_bursts[,2] == F),]
  db_alphas <- sort(db[,1])
  db_col_names <- numeric(length(db_alphas))
  for (i in 1:length(db_alphas)) {
    db_col_names[i] <- paste("db, a =", db_alphas[i])
  }
  
  jumps <- n_bursts[(n_bursts[,1] > 0) & (n_bursts[,2] == T),]
  jump_alphas <- sort(jumps[,1])
  jump_col_names <- numeric(length(jump_alphas))
  for (i in 1:length(jump_alphas)){
    jump_col_names[i] <- paste("jump, a =", jump_alphas[i])
  }
  
  colnames(output_table) <- c("Bandwidth", "Beta", "a = 0", db_col_names, jump_col_names)
  

  #### Insert values into output table
  alphas <- c(0,db_alphas,jump_alphas)
  burst_or_jump <- c(unique(table1$jump[table1$alpha == 0]), rep(F,length(db_alphas)),rep(T,length(jump_alphas)))
  
  for (a_index in 1:length(alphas)) {
     #a_index <- 4
     for (b_index in 1:length(betas)) {
       #b_index <- 1
       temp <- table1[(table1$alpha == alphas[a_index]) &
                      (table1$beta == betas[b_index]) &
                      (table1$jump == burst_or_jump[a_index]),]
       table1[(table1$alpha == db_alphas[a_index]),]
       
       if (nrow(temp)){
         for (h_index in 1:length(h_list)) {
           #h_index <- 1
           to_row <- (h_index-1)*length(betas)+b_index
           to_col <- 2+a_index
           output_table[to_row,to_col] <- round(temp[1,3+h_index],1)
         }
       }  
     }
   }
  return(output_table)
}
