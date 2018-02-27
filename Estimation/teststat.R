teststat<-function(data.mu, data.sig, ht, kern){
  # data.mu should include time / mu
  # data.sig should include time / sig
  # kern MUST contain ksq as it is being used
  
    # if whole list is passed - it will find ksq itself
  if(is.list(kern)) ksq = kern$ksq
  if(is.function(kern)) stop("kern argument should be either ksq or list containing ksq")
  
  time <- data.mu$time
  mu<-data.mu$mu
  time2 <- data.sig$time
  sig <- data.sig$sig
  
  # check if lengths makes sense
  if(length(time) != length(mu)){
    stop("length(time) != length(mu)")
  }
  if(length(time2) != length(sig)){
    stop("length(mu) != length(sig)")
  }
  
  # check if time and time2 are NOT identical
  if(!setequal(time, time2)){
    print("time.mu and time.sig differs - attempting to match")
    
    int<-intersect(time,time2)
    #adjust time, mu and sig
    mu<-mu[match(int, time)]
    sig<-sig[match(int, time2)]
    time<-int
  }
  if(length(mu) == length(sig)) print("lengths succesfully matched")
  else stop("Unable to match - clean up your input")
  
  n <- length(time)
  coef <- sqrt(ht/ksq)
  #coef <- sqrt(ht)
  
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
  res = max(data$test)
  
  return(list(start = start, end = end, tstar = res))
}

est.z_quantile<-function(mym, myrho, myalpha){
  #Implicitly uses: interpolList
  #requires function interpolListInParent()
  
  #mym: vector of 'm' values to evalute at
  #myrho: vector of 'rho' values to evaluate at
  #myalpha: chosen confidence. Needs to match possible choices in interpolList
  #interpolList: list of alphs, and list of fitted polynomials (lm-objects)
  
  #output: list of 
  #q: quantile of (Zm-am)*bm
  #qZm : raw quantile of Zm
  
  #for each pair of (mym,myrho)
  
  
  interpolListInParent()
  
  alphaIndex<-match(myalpha,interpolList$alpha)
  
  if(is.na(alphaIndex)){
    stop(paste0("Choose alpha in: ",interpolList$alpha))
  }
  
  if(length(myrho)<=1 & length(mym)<=1){
    print("Inefficient to call witth 'myrho' & 'mym' with length 1")
    mym<-rep(mym,2)
    OnedimFlag<-T
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