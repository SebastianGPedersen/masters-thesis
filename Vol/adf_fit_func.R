# adf = augmented dickey-fuller

adf.Fit.optimalP <- function(y, maxP, critF = vol.est.AICc, ...){
    # DTlist   <- vector(mode = "list", length = maxP)
    # colNameslist <- vector(mode = "list", length = maxP)
    # fitlist   <- vector(mode = "list", length = maxP)
    CritVec <- vector(mode = "numeric", length = maxP)
    
    bestCrit <- 10^100
    
    for(i in 1:maxP){
      currentFit     <- adf.Fit.SingleP(y = y, p = i)
      CritVec[i]     <- critF(currentFit, ...) 
      
      if(i%%10 == 0 ){
        print(paste0("Current p = ", i, " max p = ", maxP))
      }
      
      if(CritVec[i] + 0.01 < bestCrit ){
        bestCrit     <- CritVec[i] 
        bestFit   <- currentFit
        optimalP <- i
      }
      
    }
    return(list(optimalP = optimalP, fit = bestFit, Crit = bestCrit, CritVec = CritVec))
}

adf.Fit.SingleP <- function(y, p){

dy <- diff(y)
len <- length(dy)
strVec <- character(p)

for(i in 1:p){
  strVec[i] <- paste0("dy[-c(1:", p-i, ",",len-(i-1), ":",len,")]")
}
strVec <- gsub(pattern = "1\\:0\\,", replacement = "", strVec)
extraStr <- paste0("y[-c((1:",p,"),",len+1,")]")

fmla <- as.formula(paste0("dy[-(1:",p, ")] ~ ", paste(c(extraStr,strVec), collapse = "+")))
 
return(lm(fmla))
}


