# adf = augmented dickey-fuller

adf.Fit.optimalP <- function(y, maxP, critF = vol.est.AICc, na_Action = na.omit, ...){
    # DTlist   <- vector(mode = "list", length = maxP)
    # colNameslist <- vector(mode = "list", length = maxP)
    # fitlist   <- vector(mode = "list", length = maxP)
    CritVec <- vector(mode = "numeric", length = maxP)
    NcompleteCases <- vector(mode = "numeric", length = maxP)
    testStats <- vector(mode = "numeric", length = maxP)
    
    bestCrit <- 10^100
    
    for(i in 1:maxP){
      currentFit         <- adf.Fit.SingleP(y = y, p = i, na_Action = na_Action)
      CritVec[i]         <- critF(currentFit, ...) 
      NcompleteCases[i]  <- length(y)- length((currentFit$na.action))
      testStats[i]       <- coef(summary(currentFit))[1,1]/coef(summary(currentFit))[1,2]
      
      
      
      
      
      if(i%%10 == 0 ){
        print(paste0("Current p = ", i, " max p = ", maxP))
      }
      
      if(CritVec[i] + 0.01 < bestCrit ){
        bestCrit     <- CritVec[i] 
        bestFit   <- currentFit
        optimalP <- i
      }
      
    }
    return(list(optimalP = optimalP, fit = bestFit, Crit = bestCrit, CritVec = CritVec, NcompleteCases = NcompleteCases,  testStats = testStats))
}

adf.Fit.SingleP <- function(y, p, na_Action){

dy <- diff(y)
len <- length(dy)
strVec <- character(p)

for(i in 1:p){
  strVec[i] <- paste0("dy[-c(1:", p-i, ",",len-(i-1), ":",len,")]")
}
strVec <- gsub(pattern = "1\\:0\\,", replacement = "", strVec)
extraStr <- paste0("y[-c((1:",p,"),",len+1,")]")

fmla <- as.formula(paste0("dy[-(1:",p, ")] ~ -1 +", paste(c(extraStr,strVec), collapse = "+")))
 
return(lm(fmla, na.action = na_Action))
}


