

#Fit <- sim.BSS.Fit()
# test <- sim.BSS(1:10, 2, "Gamma")

sim.BSS.equidist_times <- function(Nsteps = 100, mat = 6.5/(24*7*52)){
  return(seq(from = 0, to = mat, length.out = Nsteps + 1)) # time in years
}

sim.BSS <- function(hVec, nPaths, S0 = 1, mu_add = 0, type = "Gamma", Fit){
  time_points <- hVec[-1] #For output
  # Assums hVec[1] = 0
  if(missing(Fit)){
    Fit <- sim.BSS.Fit()
  }
  
  hVec <- hVec * (60*24*7*52)/Fit$bvS_List$bucketLengthInMinutes # time in buckets
  dt   <- diff(hVec)
  hVec <- hVec[-1] #Remove time 0
  
  Vol <- sim.BSS.Vol(hVec, nPaths, type, Fit$alpha, Fit$memory_param, Fit$log_c_sigma, Fit$nu, Fit$bvS_List)
  dW  <- replicate(n = nPaths, rnorm(length(hVec), 0, 1), simplify = T)
  lnS <-  apply(mu_add * dt + Vol*sqrt(dt)*dW, 2, cumsum)
  # lnS <-  apply((mu_add - (Vol^2)/2) * dt + Vol*sqrt(dt)*dW, 2, cumsum)
  #S   <- S0 * exp(lnS)
  # Transpose S and Vol to get: Rows = Paths, Cols = Steps
  #S <- rbind(rep(S0, nPaths), S)
  return(list(time = time_points,X = t(lnS), vol = t(Vol)))
  
}

sim.BSS.Vol <- function(hVec, nPaths, type, alpha, memory_param, log_c_sigma, nu, bvS_List){
  
  if(type == "Power"){
    BSS_Cor <- BSS.power_cor
    BSS_Cov <- BSS.power_cov
    BSS_Var <- BSS.power_var
  } else if(type == "Gamma"){
    BSS_Cor <- BSS.gamma_cor
    BSS_Cov <- BSS.gamma_cov
    BSS_Var <- BSS.gamma_var
  } else {
    stop("Unsupported type")
  }
  
  distMat <- outer(hVec, hVec, FUN = function(x,y) abs(x-y))
  # distMat
  covMat <- distMat
  covMat[upper.tri(distMat)] <- BSS_Cov(hVec = distMat[upper.tri(distMat)], alpha = alpha, memory_param = memory_param)
  
  covMat[lower.tri(covMat)] <- t(covMat)[lower.tri(covMat)]
  # covMat # Now symmetric
  diag(covMat) <- BSS_Var(alpha = alpha, memory_param = memory_param)
  # covMat
  n <- length(hVec)
  covMat <- covMat * nu^2
  genMat <- t(chol(covMat)) 
  
  #Save covariance matrix
  setwd(Sys.getenv("masters-thesis"))
  save(genMat, file = "../Personal/CovMat.Rda")
  
  log_Deseason_Sigma_Paths <- replicate(n = nPaths, genMat %*% rnorm(n, 0, 1), simplify = T)
  seasonal_Component <- sim.BSS.Seasonality(hVec, bvS_List)
  
  return(exp(log_c_sigma + log_Deseason_Sigma_Paths + seasonal_Component))
}

save.BSS.cov <- function(hVec, nPaths, S0 = 1, mu_add = 0, type = "Gamma", Fit){

  # Assums hVec[1] = 0
  if(missing(Fit)){
    Fit <- sim.BSS.Fit()
  }
  
  hVec <- hVec * (60*24*7*52)/Fit$bvS_List$bucketLengthInMinutes # time in buckets
  dt   <- diff(hVec)
  hVec <- hVec[-1] #Remove time 0  
  
  ###Change to sim.BSS.Vol types
  alpha <- Fit$alpha
  memory_param <- Fit$memory_param
  log_c_sigma <- Fit$log_c_sigma
  nu <- Fit$nu
  bvS_List <- Fit$bvS_List
  
  ###Run Sim.BSS.vol both without copy of covariance matrix
  if(type == "Power"){
    BSS_Cor <- BSS.power_cor
    BSS_Cov <- BSS.power_cov
    BSS_Var <- BSS.power_var
  } else if(type == "Gamma"){
    BSS_Cor <- BSS.gamma_cor
    BSS_Cov <- BSS.gamma_cov
    BSS_Var <- BSS.gamma_var
  } else {
    stop("Unsupported type")
  }
  
  covMat <- outer(hVec, hVec, FUN = function(x,y) abs(x-y))

  covMat[upper.tri(covMat)] <- BSS_Cov(hVec = covMat[upper.tri(covMat)], alpha = alpha, memory_param = memory_param)
  
  covMat[lower.tri(covMat)] <- t(covMat)[lower.tri(covMat)]
  # covMat # Now symmetric
  diag(covMat) <- BSS_Var(alpha = alpha, memory_param = memory_param)
  # covMat
  n <- length(hVec)
  covMat <- covMat * nu^2
  covMat <- t(chol(covMat)) 
  
  #Save covariance matrix

  setwd(Sys.getenv("masters-thesis-data"))
  save(covMat, file = "BSS_sim/covMat.Rda")

}



sim.BSS.Seasonality <- function(hVec, bvS_List) {
  newDT <- data.table(PredictTimes = hVec)
  FFF.initFFF(DT = newDT, P = bvS_List$optimalP, IntradayCol = "PredictTimes", Period = bvS_List$Period) #Edits by reference
  return(predict(bvS_List$FFF$fit, newDT))
  
  # plot(FFF.singlePeriodFit(res$fit, bvSDTfff, "IntraDayBucket"), type = "l")
  # PredictTimes <- seq(0.1,78, by = 0.1)
  # newDT <- data.table(PredictTimes = PredictTimes)
  # FFF.initFFF(DT = newDT, P = res$optimalP, IntradayCol = "PredictTimes", Period = max(bvSDTfff$IntraDayBucket))
  # newDT
  # pred <- predict(res$fit, newDT)
  # lines(PredictTimes,  pred, col = "red")
}

sim.BSS.Fit <- function(bucketLengthInMinutes = 5, type = "Gamma", maxLag = 125, m = 6, dt, bvS_List) {
  ########## DATA
  if(missing(dt)){
    dt <- BV.get_SPY_data()
  }
  ########### Vol estimation
  if(missing(bvS_List)){
    bvS_List <- BV.data_deseason_BV_Func(dt = dt, bucketLengthInMinutes = bucketLengthInMinutes)
  }
  ##### BSS STEP 1
  #c_ for current
  c_varDT <- PR.variogram_prep_DT(bvSDTfff = bvS_List$bvSDTfff)
  c_alpha <- PR.est.alpha(variogramDT = c_varDT, m = m, bucketLengthInMinutes = bucketLengthInMinutes, OLS = T)
  c_alpha <- unname(c_alpha)
  
  
  ##### BSS STEP 2
  if(type == "Power"){
    BSS_Cor <- BSS.power_cor
    BSS_Cov <- BSS.power_cov
    BSS_Var <- BSS.power_var
  } else if(type == "Gamma"){
    BSS_Cor <- BSS.gamma_cor
    BSS_Cov <- BSS.gamma_cov
    BSS_Var <- BSS.gamma_var
  } else {
    stop("Unsupported type")
  }
  
  # if(SPY_Bool)
    emp <- acf(x = bvS_List$bvSDT$LogVolnCorrect, lag.max = maxLag, plot = F)$acf[,,1][-1]
  # else { # Have not made else
  #   stop("Bitcoin not supported yet")
  #   #emp <- PR.acfFunc(bvSDT[(1 + i-1):s[i]], lags = 1:c_m_max, logVolCol = "LogVolnCorrect")$acf
  # }
  
  #only until negative
  if (min(emp) < 0) {
    # ERR check
    warning("max lag used:")
    #change
    maxLag <- match(1, emp < 0) - 1 # get last index before negative
    if(is.na(maxLag)){
      stop("All acf negative")
      break
    }
    # if(SPY_Bool)
    emp <- acf(x = bvS_List$bvSDT$LogVolnCorrect, lag.max = maxLag, plot = F)$acf[,,1][-1]
    # else { # Have not made else
    #   stop("Bitcoin not supported yet")
    #   #emp <- PR.acfFunc(bvSDT[(1 + i-1):s[i]], lags = 1:c_m_max, logVolCol = "LogVolnCorrect")$acf
    # }
  }
    
    # if(SPY_Bool){
    h <- 1:maxLag
    # } else {

    # }
    
  if(type == "Power"){
    limit <- 0.501
  } else {
    limit <- 0.000001
  }
  if(type == "Power"){
    optimRes <- optim(c(0.65), function(x, hVec, emp, BSS_Cor){optimF(c(c_alpha, x), hVec, emp, BSS_Cor)}, hVec = h, emp = emp, BSS_Cor=BSS_Cor , method = "L-BFGS-B",  lower = c(limit), upper = c(2000), control = list(maxit = 200, factr = 1e13))
  } else {
    # if(SPY_Bool) {
    optimRes <- optim(c(0.5), function(x, hVec, emp, BSS_Cor){optimF(c(c_alpha, x), hVec, emp, BSS_Cor)}, hVec = h, emp = emp, BSS_Cor=BSS_Cor , method = "BFGS")  
  }
    
    c_memory_param <- optimRes$par 
    
    # ERR check
    if(optimRes$convergence != 0){
      stop("No convergence in memory parameter estimation")
    }
    
    if(abs(c_alpha)>1/2 || c_memory_param < limit){
      stop("parameter space broken")
      if(!ParameterSpaceBroken){
        First_ps_offense <- i
        ParameterSpaceBroken <- T
      }
    }
    
    # BSS STEP 3
    log_c_sigma <- mean(bvS_List$bvSDT$LogVolnCorrect)
    xi <- bvS_List$bvSDT$LogVolnCorrect-log_c_sigma # Only temporary parameter
    
    # BSS STEP 4 
    secondMoment <- mean(xi^2) #equal to variance, as 0 first moment
    nu <- sqrt(secondMoment/BSS_Var(alpha = c_alpha, memory_param = c_memory_param))
    
    return(list(bvS_List = bvS_List, alpha = c_alpha, memory_param = c_memory_param, log_c_sigma = log_c_sigma, nu = nu))
    
}
