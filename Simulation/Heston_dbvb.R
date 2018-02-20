source(paste(Sys.getenv("masters-thesis"),"Simulation/Heston.R",sep="/"))


sim.all <- function(settings,burst_time,interval_length,alpha,beta,c_1,c_2) {
  
  #Initializations
  time = 1:settings$steps
  
  X = matrix(nrow = settings$Npaths, ncol = settings$steps)
  Y = matrix(nrow = settings$Npaths, ncol = settings$steps)
  vol = matrix(nrow = settings$N, ncol = settings$steps)
  
  #Define burst interval
  burst_begin = burst_time-interval_length/2
  burst_end = burst_time+interval_length/2
  mat = settings$mat #save maturity
  steps = settings$Nsteps #save nsteps
  
  if (burst_begin < 0 | burst_end > ttm) {
    stop("Burst interval is not contained in simulation interval")
  }
  
  
  # ----------------------- Heston simulation before burst interval ----------------------
  
    #Arguments for sim.heston (take steps until inside interval)
    Nsteps_b = ceil((burst_begin/mat)*steps)
    mat_b = mat*(Nsteps_b/steps)
    
    settings$Nsteps = Nsteps_b
    settings$mat = mat_b
    
    #Points before interval
    [time_b,Y_b,X_b,vol_b] <- sim.heston(settings)
  
  
  # ----------------------- Heston, VB and DB+VB inside interval ----------------------
  
    #Arguments for sim functions
    burst_time = burst_time-settings$mat #Burst time from beginning of inside interval
    Nsteps_interval = floor(burst_end/mat*steps+1) - Nsteps_b #Because of "]", we take another step inside on boundary
    mat_interval = mat*(Nsteps_interval/steps)

    settings$Nsteps = Nsteps_interval
    settings$mat = mat_interval    
    
    #Initialize vol and X from before interval
    startvol = vol_b[,ncol(vol_b)]
    X_init = X_b[,ncol(X_b)]
    
    #Calculate common rnorms
    dW = replicate(steps-1,rnorm(N,0,1))
    epsilon = rep(steps,rnorm(N,0,1))
    
    args = list(dW=dW,X_init=X_init,epsilon=epsilon,startvol=startvol)
    
    #Points inside interval
    list(time_H,Y_H,X_H,vol_H) = sim.heston(settings,arg)
    
    list(time_vb,Y_vb,X_vb,vol_vb) = 
      sim.vb(burst_time = burst_time,Nsteps = settings$Nsteps,mat = settings$mat,
             startvol = startvol,X_init = X_init,dW = dW,epsilon = epsilon,
             c_2 = c_2,beta = beta,gamma = gamma)
    
    list(time_dbvb,Y_dbvb,X_dbvb,vol_dbvb) = 
      sim.dbvb(burst_time = burst_time,Nsteps = settings$Nsteps,mat = settings$mat,
               startvol = startvol,X_init = X_init,dW = dW,epsilon = epsilon,
               c_2 = c_2,beta = beta,gamma = gamma, c_1=c_1, alpha=alpha)

  
  # ----------------------- Heston again outside interval ----------------------
    
    #Arguments for sim.heston
    Nsteps_a = steps-Nsteps_b-Nsteps_interval
    mat_a = mat-mat_b-mat_interval
    
    settings$Nsteps = Nsteps_a
    settings$mat = mat_a
    
    #Initialize vol and X from before interval
    startvol = vol_b[,ncol(vol_b)]
    X_init = X_b[,ncol(X_b)]
    
    #Calculate common rnorms
    dW = replicate(steps-1,rnorm(N,0,1))
    epsilon = rep(steps,rnorm(N,0,1))
    
    args_H = list(dW=dW,epsilon=epsilon,startvol=vol_H[,ncol(vol_H)],X_init = X_H[,ncol(X_H)])
    args_db = list(dW=dW,epsilon=epsilon,startvol=vol_vb[,ncol(vol_vb)],X_init = X_vb[,ncol(X_vb)])
    args_dbvb = list(dW=dW,epsilon=epsilon,startvol=vol_dbvb[,ncol(vol_dbvb)],X_init = X_dbvb[,ncol(X_dbvb)])
    
    #Points after interval
    list(time_H_a,Y_H_a,X_H_a,vol_H_a) = sim.heston(settings,args_H)
    list(time_vb_a,Y_vb_a,X_vb_a,vol_vb_a) = sim.heston(settings,args_db)
    list(time_dbvb_a,Y_dbvb_a,X_dbvb_a,vol_dbvb_a) = sim.heston(settings,args_dbvb) 


    # ----------------------- Wrap together and return ----------------------
    time = 0:steps
    
    Heston = 
    vb = 
    dbvb =
    
    return(list(Heston = Heston,vb=vb,dbvb=dbvb))
}





sim.vb <- function(burst_time,Nsteps,mat,startvol,X_init,dW,epsilon,c_2,beta,gamma){
  #Also returns values at time zero
  
  time = 0:Nsteps
  X = matrix(nrow = N, ncol = Nsteps+1)
  Y = matrix(nrow = N, ncol = Nsteps+1)
  vol = matrix(nrow = N, ncol = Nsteps+1)
  dt = mat/Nsteps
  
  
  #init
  vol[,1] = startvol
  X[,1] = X_init
  Y[,1] = X_init #Y_init doesn't matter, we remove it anyway
  
  
  for(i in 2:(Nsteps+1)){
    vol[,i] = c_2*sign(dt*(i-1)-burst_time)/abs(dt*(i-1)-burst_time)^beta #calulate vol
    NS = dW[,i-1]
    
    X[,i] = sqrt(vol[,i])*NS
    
    omega = gamma*sqrt(vol[,i])/sqrt(Nsteps)     # n corresponds to steps and not repetitions N? #should vol be i-1? No?
    Y[,i] = X[,i] + omega * epsilon[,i]
  }
  return(list(time = time, Y = Y, X = X, vol = vol))
}


sim.dbvb <- function(burst_time,Nsteps,mat,startvol,X_init,dW,epsilon,c_2,beta,gamma,c_1,alpha){
  #Also returns values at time zero
  
  time = 0:Nsteps
  X = matrix(nrow = N, ncol = Nsteps+1)
  Y = matrix(nrow = N, ncol = Nsteps+1)
  vol = matrix(nrow = N, ncol = Nsteps+1)
  dt = mat/Nsteps
  
  
  #init
  vol[,1] = startvol
  X[,1] = X_init
  Y[,1] = X_init #Y_init doesn't matter, we remove it anyway
  
  
  for(i in 2:(Nsteps+1)){
    drift = c_1*sign(dt*(i-1)-burst_time)/abs(dt*(i-1)-burst_time)^alpha #calculate drift
    vol[,i] = c_2*sign(dt*(i-1)-burst_time)/abs(dt*(i-1)-burst_time)^beta #calulate vol
    NS = dW[,i-1]
    
    X[,i] = drift*dt+sqrt(vol[,i])*NS
    
    omega = gamma*sqrt(vol[,i])/sqrt(Nsteps)     # n corresponds to steps and not repetitions N? #should vol be i-1? No?
    Y[,i] = X[,i] + omega * epsilon[,i]
  }
  return(list(time = time, Y = Y, X = X, vol = vol))
}