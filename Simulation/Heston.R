
# IT'S NOT WHO YOU ARE UNDERNEATH BUT WHAT YOU DO THAT DEFINES YOU


# List of needed input that can be sent to simulation function
sim.setup <- function(kappa=5, theta=0.0225, xi = 0.4, rho = -0.5, omega = 1.6*10^(-5),
                   mat = 6.5/(24*7*52), Nsteps = 23400, Npath = 1000){
  list(kappa = kappa, theta = theta, xi = xi, rho = rho, omega = omega, mat = mat, Nsteps = Nsteps, Npath = Npath)
}

# Heston simulation (no scheduled bursts)
sim.heston<-function(settings){
  
  N = settings$Npath
  mat = settings$mat
  steps = settings$Nsteps
  kappa = settings$kappa
  theta = settings$theta
  xi = settings$xi       ##vol of vol
  rho = settings$rho
  omega = settings$omega
  
  dt = mat/steps #dt is in years
  time = 0:steps*dt # maybe time should be relative (0-1)
  #time = 0:(steps)/(steps)
  
  X = matrix(nrow = N, ncol = steps+1)
  Y = matrix(nrow = N, ncol = steps+1)
  vol = matrix(nrow = N, ncol = steps+1) #matrix if we want to save values along the way
  call = numeric(N) # testing purpose
  
  X[, 1] = 0
  vol[, 1] = rgamma(N, 2*kappa*theta/xi^2, 2*kappa/xi^2)
  #Y[, 1] = X[,1] + gamma*sqrt(vol[,1])/sqrt(steps)*rnorm(N,0,1) #Changed from vol to sqrt(vol) /Seb 20.02.18
  Y[ ,1] = X[, 1] + omega*rnorm(N,0,1)
  
  for(i in 2:(steps+1)){
    NS = rnorm(N,0,1)
    NV = rho*NS + sqrt(1-rho^2)*rnorm(N,0,1) #From StatÃ (Olivier) Theorem I.5 or Graphical example 1.20 
    
    #X[,i] =   X[,i-1] -1/2*vol[,i-1]*dt+ sqrt(vol[,i-1])*sqrt(dt)*NS        # real heston
    #X[,i] =   X[,i-1] + X[,i-1]*sqrt(vol[,i-1])*sqrt(dt)*NS                 # non-log
    X[,i] =   X[,i-1] + sqrt(vol[,i-1])*sqrt(dt)*NS                          # logs (driftless)
    
    x = vol[,i-1] + kappa*(theta - vol[,i-1])*dt
    y = sqrt(  log(  (xi^2*vol[,i-1]*dt)/(x^2) + 1  )  )
    
    vol[, i] = x*exp(-0.5*y^2+y*NV  )
    
    #Observed Y
    #omega = gamma*sqrt(vol[,i])/sqrt(steps)     # n corresponds to steps and not repetitions N? #should vol be i-1? No?
    Y[,i] = X[,i] + omega * rnorm(N,0,1)
    
  }
  return(list(time = time, Y = Y, X = X, vol = vol))
}

# Helper functions to streamline simulation object to estimation object setup
# (list of matrices from sim -> vector of lists of vectors)

sim.path <- function(path, sim.data, compact=T){
  # Extracts the data list of a single path from the simulation data

  m <- dim(sim.data$Y)[1] # Npath
  if(path > m) stop("out of bounds - max is: ", m)
  if(compact){
    return(list(time = sim.data$time, Y = sim.data$Y[path,]))
  }
  else{
    return(list(time = sim.data$time, Y = sim.data$Y[path,], X = sim.data$X[path,], vol = sim.data$vol[path,]))
  }
}

# Heston simulation (no scheduled bursts)
sim.heston.uneven<-function(settings, seed = 0){
  # simulates similarily to the real-life distribution of data points
  N = settings$Npath
  mat = settings$mat
  steps = settings$Nsteps
  kappa = settings$kappa
  theta = settings$theta
  xi = settings$xi       ##vol of vol
  rho = settings$rho
  omega = settings$omega
  
  # BEGIN TIME ACROBATIQUE
  trades_dist<-readRDS("Simulation/trades_dist.Rda")
  #normalize
  trades_dist_norm<-trades_dist/(sum(trades_dist))
  
  # find number of trades per 'bucket'
  trades<-floor(trades_dist_norm*steps)
  trades[length(trades)] <- trades[length(trades)] + max(steps - sum(floor(trades_dist_norm*steps)),0)
  
  # Simulate dist within each bucket
  time <- runif(trades[1])
  for(i in 2:length(trades)){
    time<-c(time,runif(trades[i])+(i-1))
  }
  time<-c(0,mat*sort(time)/length(trades))
  # END TIME ACROBATIQUE
  
  dt <- diff(time) # should be length(step)
  
  if(!seed == 0){
    set.seed(seed)
  }
  
  X = matrix(nrow = N, ncol = steps+1)
  Y = matrix(nrow = N, ncol = steps+1)
  vol = matrix(nrow = N, ncol = steps+1)
  
  X[, 1] = 0
  vol[, 1] = rgamma(N, 2*kappa*theta/xi^2, 2*kappa/xi^2)
  Y[ ,1] = X[, 1] + omega*rnorm(N,0,1)
  
  for(i in 2:(steps+1)){
    NS = rnorm(N,0,1)
    NV = rho*NS + sqrt(1-rho^2)*rnorm(N,0,1)
    
    X[,i] =   X[,i-1] + sqrt(vol[,i-1])*sqrt(dt[i-1])*NS
    
    x = vol[,i-1] + kappa*(theta - vol[,i-1])*dt[i-1]
    y = sqrt(  log(  (xi^2*vol[,i-1]*dt[i-1])/(x^2) + 1  )  )
    
    vol[, i] = x*exp(-0.5*y^2+y*NV  )
    
    #Observed Y
    Y[,i] = X[,i] + omega * rnorm(N,0,1)
    
  }
  return(list(time = time, Y = Y, X = X, vol = vol))
}
