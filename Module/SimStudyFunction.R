study<-function(setting, hd, hv, t.index, conf, burstsetting){
  # The study simulates a heston, adds vol burst and evaluates the critical T* value
  # to replicate table 1
  
  alpha <- burstsetting$alpha
  beta <- burstsetting$beta
  burst_time <- burstsetting$burst_time
  interval_length <- burstsetting$interval_length
  c_1 <- burstsetting$c_1
  c_2 <- burstsetting$c_2
  
  # simulation
  heston<-sim.heston(setting)
  
  # pathwise
  N<-dim(sim$Y)[1]
  Tstar <- numeric(N)
  rhom <- numeric(N)
  rhorho<-numeric(N)
  
  Tstarraw <- numeric(N)
  rhomraw <- numeric(N)
  rhorhoraw<-numeric(N)
  
  Tstarvb <- numeric(N)
  rhomvb <- numeric(N)
  rhorhovb<-numeric(N)
  
  Tstardb <- numeric(N)
  rhomdb <- numeric(N)
  rhorhodb<-numeric(N)
  for(mode in 1:3){
    if(mode == 2){
      sim <- sim.addvb(sim,    burst_time = burst_time, interval_length = interval_length,
                       c_2 = c_2, beta  = beta)
    }
    if(mode == 3){
      sim <- sim.adddb(sim, burst_time = burst_time, interval_length = interval_length,
                       c_1 = c_1,  alpha = alpha)
    }
    for(i in 1:N){
      # Extract
      simpath<-sim.path(i, sim)
      
      #data<-everyOther(simpath) lav det 
      
      # Estimation of mu/sig
      mu<-est.mu(simpath, hd, kern.leftexp, t.index = t.index)
      sig <- est.sigma(simpath, hv, kern.leftexp, kern.parzen, t.index = t.index, lag = "auto")
      
      # Calculate T
      Tstat<-teststat(mu, sig, hd, hv)
      
      # Calculate T*
      Tstar[i]<-tstar(Tstat)$tstar
      
      # fit rho
      rho <- est.rho(Tstat$test)
      rhom[i] <- rho$m
      rhorho[i] <- rho$rho
      
    }
    if(mode == 1){
      Tstatraw <- Tstat
      rhomraw <- rhom
      rhorhoraw <- rhorho
    }
    else if(mode == 2){
      Tstatvb <- Tstat
      rhomvb <- rhom
      rhorhovb <- rhorho
    }
    else{
      Tstatdb <- Tstat
      rhomdb <- rhom
      rhorhodb <- rhorho
    }
  }
  # plug in rest (q-quantile)
  zraw<-est.z_quantile(rhomraw, rhorhoraw, conf)$qZm
  zvb<-est.z_quantile(rhomvb, rhorhovb, conf)$qZm
  zdb<-est.z_quantile(rhomdb, rhorhodb, conf)$qZm
  resraw<-mean(Tstarraw>=zraw)
  resvb<-mean(Tstarvb>=zvb)
  resdb<-mean(Tstardb>=zdb)
  
  input <- list(setting = setting, hd = hd, hv = hv, t.index = t.index,
                conf = conf, burstsetting   = burstsetting)
  
  output <- list(raw = resraw, vb = resvb, db = resdb,
                 Tstarraw = Tstarraw, Tstarvb = Tstarvb, Tstardb = Tstardb)
  
  return(list(input = input, output = output))
}