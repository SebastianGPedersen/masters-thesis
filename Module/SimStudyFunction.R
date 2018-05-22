study<-function(setting, burstsetting, hd, hv, t.index, conf){
  # The study simulates a heston, adds vol burst and evaluates the critical T* value
  # to replicate table 1
  
  alpha <- burstsetting$alpha;  beta <- burstsetting$beta;
  burst_time <- burstsetting$burst_time;  interval_length <- burstsetting$interval_length;
  c_1 <- burstsetting$c_1;  c_2 <- burstsetting$c_2
  
  # simulation
  heston<-sim.heston(setting)
  
  # pathwise
  N<-dim(heston$Y)[1]
  
  Tstar <- numeric(N);  rhom <- numeric(N);  rhorho<-numeric(N);
  
  Tstarraw   <- numeric(N);  rhomraw   <- numeric(N);  rhorhoraw   <- numeric(N);
  Tstarvb    <- numeric(N);  rhomvb    <- numeric(N);  rhorhovb    <- numeric(N);
  Tstardb    <- numeric(N);  rhomdb    <- numeric(N);  rhorhodb    <- numeric(N);
  TstarJ     <- numeric(N);  rhomJ     <- numeric(N);  rhorhoJ     <- numeric(N);
  TstarJvb   <- numeric(N);  rhomJvb   <- numeric(N);  rhorhoJvb   <- numeric(N);
  
  for(mode in 1:5){
    if(mode == 1){
      sim <- heston
    }
    if(mode == 2){
      sim <- sim.addvb(heston,    burst_time = burst_time, interval_length = interval_length,
                       c_2 = c_2, beta  = beta)
    }
    if(mode == 3){
      sim <- sim.addvb(heston,    burst_time = burst_time, interval_length = interval_length,
                       c_2 = c_2, beta  = beta)
      sim <- sim.adddb(sim, burst_time = burst_time, interval_length = interval_length,
                       c_1 = c_1,  alpha = alpha)
    }
    if(mode == 4){
      sim <- sim.addjump(heston, burst_time = burst_time, interval_length = interval_length,
                         c_1 = c_1, alpha = alpha)
    }
    if(mode == 5){
      sim <- sim.addvb(heston,    burst_time = burst_time, interval_length = interval_length,
                       c_2 = c_2, beta  = beta)
      sim <- sim.addjump(sim, burst_time = burst_time, interval_length = interval_length,
                         c_1 = c_1, alpha = alpha)
    }
    if(mode > 5){
      stop("Error in mode")
    }
    #N is path
    data <- list(time = sim$time, Y = t(diff(t(sim$Y))) )
    # Estimation of mu/sig
    mu  <-est.mu.mat.next(data, hd = hd, t.index = t.index)
    sig <- est.sigma.mat.next(data, hv = hv, t.index = t.index, lag = 10)
    
    # Calculate T
    Tstat<-teststat(mu, sig, hd, hv) #matrix output
    
    # Calculate T*
    for(path in 1:N){
      Tstar[path] <- max(abs(Tstat$test[path,]))
      # fit rho
      rho <- est.rho(Tstat$test[path,])
      rhom[path] <- rho$m
      rhorho[path] <- rho$rho
    }
    
    # Assign to correct output
    if(mode == 1){
      Tstarraw <- Tstar
      #Tstatraw <- Tstat
      rhomraw <- rhom
      rhorhoraw <- rhorho
    }
    else if(mode == 2){
      Tstarvb <- Tstar
      #Tstatvb <- Tstat
      rhomvb <- rhom
      rhorhovb <- rhorho
    }
    else if(mode == 3){
      Tstardb <- Tstar
      #Tstatdb <- Tstat
      rhomdb <- rhom
      rhorhodb <- rhorho
    }
    else if(mode == 4){
      TstarJ <- Tstar
      #TstatJ <- Tstat
      rhomJ <- rhom
      rhorhoJ <- rhorho
    }
    else{
      TstarJvb <- Tstar
      #TstatJvb <- Tstat
      rhomJvb <- rhom
      rhorhoJvb <- rhorho
    }
  }
  # plug in rest (q-quantile)
  zraw<-est.z_quantile(rhomraw, rhorhoraw, conf)$qZm #raw
  zvb<-est.z_quantile(rhomvb, rhorhovb, conf)$qZm    #vb
  zdb<-est.z_quantile(rhomdb, rhorhodb, conf)$qZm    #db
  zJ<-est.z_quantile(rhomJ, rhorhoJ, conf)$qZm    #J
  zJvb<-est.z_quantile(rhomJvb, rhorhoJvb, conf)$qZm    #vbJ
  resraw<-mean(Tstarraw>=zraw)
  resvb<-mean(Tstarvb>=zvb)
  resdb<-mean(Tstardb>=zdb)
  resJ<-mean(TstarJ>=zJ)
  resJvb<-mean(TstarJvb>=zJvb)
  
  input <- list(setting = setting, hd = hd, hv = hv, t.index = t.index,
                conf = conf, burstsetting   = burstsetting)
  
  output <- list(raw = resraw, vb = resvb, db = resdb,
                 J = resJ, Jvb = resJvb,
                 Tstarraw = Tstarraw, Tstarvb = Tstarvb, Tstardb = Tstardb,
                 TstarJ = TstarJ, TstarJvb = TstarJvb)
  
  return(list(input = input, output = output, debug = heston))
}