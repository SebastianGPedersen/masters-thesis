est.PreAverage<-function(logRet_Data, K){
  if((!(K%%2==0))&(K>1)){
    stop("Needs K even")
  }
  
  
  N <- length(logRet_Data) 
  lOut <- N-K+1
  r <- rep(NA, lOut)
  
  for(i in 1:(lOut)){
    
    r[i] <- sum(logRet_Data[i+(K/2):(K-1)])-sum(logRet_Data[i+(0):(K/2-1)])
  }  
  
  return(r/K)
}
