FFF.estOptimalP<- function(DT, maxP, dataCol,  IntradayCol, dayCol, dailyInteraction = F, critF = FFF.AICc, dailyParameter = T, minP = 1, ...){
  # DTlist   <- vector(mode = "list", length = maxP)
  # colNameslist <- vector(mode = "list", length = maxP)
  # fitlist   <- vector(mode = "list", length = maxP)
  CritVec <- vector(mode = "numeric", length = maxP)
  

  bestCrit <- 10^100
  
  for(i in minP:maxP){
    currentDT           <- copy(DT)
    currentColNames     <- FFF.initFFF(currentDT, i, IntradayCol) #changes DT by reference
    
    if(dailyInteraction){
      fmla         <- as.formula(paste(dataCol , "~ 1 +", dayCol, "*(", paste(currentColNames,collapse = "+"), ")"))
    } else {
      if(dailyParameter){
        fmla       <- as.formula(paste(dataCol , "~ 1 +", dayCol, "+", paste(currentColNames,collapse = "+")))
      } else {
        fmla       <- as.formula(paste(dataCol , "~ 1 +", paste(currentColNames,collapse = "+")))
      }
    }

    currentFit   <- lm(fmla, data = currentDT)
    CritVec[i]     <- critF(currentFit, ...) 
    
    if(i%%10 == 0 ){
      print(paste0("Current p = ", i, " max p = ", maxP))
    }
    
    if(CritVec[i] + 0.01 < bestCrit ){
      bestCrit     <- CritVec[i] 
      bestDT       <- currentDT
      bestFit   <- currentFit
      bestColNames <- currentColNames
      optimalP <- i
    }
    
  }
  return(list(optimalP = optimalP, DT = bestDT, fit = bestFit, ColNames = bestColNames, Crit = bestCrit, CritVec = CritVec))
}

FFF.initFFF <- function(DT, P, IntradayCol){
  #EDITS DT BY REFERENCE
  #data.table DT with cols
    #IntradayCol (column name as string): Col with repeating numbers 1:x, with x denoting a number of intervals the day has been split into (5min buckets, 6.5h day = 78)
  #P number of sin & cos terms to include
  
  # vol <- DT$volCol #extract only once
  intraDayVec <- DT[[IntradayCol]]
  Period <- max(intraDayVec)
  
  sinTerm <- function(p){
    #NOT vectorized in p
    sin(intraDayVec*2*pi*p/Period)
  }
  
  cosTerm <- function(p){
    #NOT vectorized in p
    cos(intraDayVec*2*pi*p/Period)
  }
  sinL <- vector(mode = "list", length = P)
  cosL <- vector(mode = "list", length = P)
  
  for(i in 1:P){
    sinL[[i]] <- sinTerm(i)
    cosL[[i]] <- cosTerm(i)
  }
  
  polyTermL <- vector(mode = "list", length = 2)
  polyTermL[[1]] <- intraDayVec
  polyTermL[[2]] <- intraDayVec^2
  newColNames <- c("poly1", "poly2", paste0("sin", 1:P),paste0("cos", 1:P))
  
  #Adds new cols
  DT[, (newColNames):=c(polyTermL,sinL,cosL)]
  
  #for easy construction of formula
  return(newColNames)
}

FFF.singlePeriodFit <- function(fit, DT, IntradayCol){
  #### VERY specific requirements - use with caution.

  #fit = fitted FFF
  #DT = same DT as used for fit
  #newColNames = output from FFF.initFFF
  fit$fitted.values[1:max(DT[[IntradayCol]])]#-fit$coefficients["(Intercept)"]
  # (as.matrix(fit$model[newColNames])%*%as.matrix(fit$coefficients[newColNames]))[1:max(DT[[IntradayCol]])]
}


FFF.AICc<-function(fit){
  logL <- logLik(fit)
  L <- logL[1]
  k <- attributes(logL)$df
  n <- attributes(logL)$nobs
  return(2*(k-L+(k^2+k)/(n-k-1)))
}

FFF.AIC<-function(fit){
  logL <- logLik(fit)
  L <- logL[1]
  k <- attributes(logL)$df
  n <- attributes(logL)$nobs
  return(2*(k-L))
}