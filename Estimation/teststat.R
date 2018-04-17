teststat<-function(data.mu, data.sig, hd, hv, kern=kern.leftexp){
  # data.mu should include time / mu
  # data.sig should include time / sig
  
  time <- data.mu$time
  mu<-data.mu$mu
  sig <- data.sig$sig
  
  coef <- sqrt(hd/kern$ksq)    #sqrt(hd)*sqrt(hv) <- our previous version
  
  # calculates t-stat
  t <- coef*(mu/sqrt(sig))
  
  # returns list of the grouped and the raw
  return(list(time = time, test = t))
}

# change this
tstar<-function(data){
  # data should contain time and test
  # Need vector/list of start and end for each period (or n. of observations in)
  start = data$time[1]
  end = data$time[length(data$time)]
  res = max(abs(data$test))
  
  return(list(start = start, end = end, tstar = res))
}

test.db<-function(data, hd, hv, kern = kern.leftexp, noisefun, theta, kn){
  # data should include time | Y (log returns) | raw (before preaverage)
  # Calculates db test stat by computing mu/sig
  
  mu<-est.mu.new(data = data, hd = hd, t.index = tind, kn = k)
  sig <- est.sigma.new(data, hv=hv, t.index = tind, noisefun = est.noise.iid, theta = theta, kn = k)
  
  # Calculate T
  Tstat<-teststat(mu, sig, hd, hv, kern)
  
  return(Tstat)
}

est.z_quantile<-function(mym, myrho, myalpha){
  #Vectorized in mym and myrho (should be same length)
  
  #Implicitly uses: interpolList
  #requires function est.interpolListInParent()
  
  #mym: vector of 'm' values to evalute at
  #myrho: vector of 'rho' values to evaluate at
  #myalpha: chosen confidence. Needs to match possible choices in interpolList
  #interpolList: list of alphs, and list of fitted polynomials (lm-objects)
  
  #output: list of 
  #q: quantile of (Zm-am)*bm
  #qZm : raw quantile of Zm
  
  #for each pair of (mym,myrho)
  
  
  est.interpolListInParent()
  
  alphaIndex<-match(myalpha,interpolList$alpha)
  
  if(is.na(alphaIndex)){
    stop(paste0("Choose alpha in: ",paste0(interpolList$alpha, collapse = ", ")))
  }
  
  if( myrho < 0 || 0.999 < myrho ){
    print("Extrapolating on rho")
  }
  
  if( mym < 100 || 10000 < mym ){
    print("Extrapolating on m")
  }
  
  if(length(myrho)<=1 & length(mym)<=1){
    print("Inefficient to call with 'myrho' & 'mym' with length 1")
    mym<-rep(mym,2)
    OnedimFlag<-T
  } else {
    OnedimFlag<-F
  }
  fit<-interpolList$fittedPoly[[alphaIndex]]
  newdata<-data.frame(logm = log(mym), logrho = log(1-myrho))
  q<-predict(fit,newdata)
  #q
  
  am <- sqrt(2*log(mym));   
  bm <- am-0.5*log(pi*log(mym))/am;
  qZm <- q/am+bm;  
  #qZm
  
  if(OnedimFlag){
    q <- q[1]
    qZm <- qZm[1]
  }
  
  return(list(q=q, qZm=qZm))
}


est.interpolListInParent<-function(){
  #creates interpolList in the global environment, if it does not exsists.
  if(!exists("interPolList")){
    cd<-getwd()
    setwd(Sys.getenv("masters-thesis-data"))
    assign("interpolList", readRDS("interpolList.rds"), envir = .GlobalEnv)
    setwd(cd)
  }
}
