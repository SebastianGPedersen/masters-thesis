
# IT'S NOT WHO YOU ARE UNDERNEATH BUT WHAT YOU DO THAT DEFINES YOU


# List of needed input that can be sent to simulation function
sim.setup <- function(kappa=5, theta=0.0225, xi = 0.4, rho = -0.5, gamma = 0.5,
                   mat = 1, Nsteps = 1000, Npath = 10000){
  list(kappa = kappa, theta = theta, xi = xi, rho = rho, gamma = gamma, mat = mat, Nsteps = Nsteps, Npath = Npath)
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
  gamma = settings$gamma
  
  dt = mat/steps
  time = 1:steps
  
  X = matrix(nrow = N, ncol = steps)
  Y = matrix(nrow = N, ncol = steps)
  vol = matrix(nrow = N, ncol = steps) #matrix if we want to save values along the way
  
  startvol = rgamma(N, 2*kappa*theta/xi^2, 2*kappa/xi^2)
  
  X[, 1] = 1
  vol[, 1] = startvol
  Y[, 1] = X[,1] + gamma*vol[,1]/sqrt(steps)*rnorm(N, 0, 1) #n corresponds to steps and not repetitions N?
  
  for(i in 2:steps){
    NS = rnorm(N, 0, 1)
    NV = rho*NS + sqrt(1-rho^2)*rnorm(N,0,1)
    
    #X[,i] =   X[,i-1] + X[,i-1]*sqrt(vol[,i-1])*sqrt(dt)*NS         #non-ln x's
    X[,i] =   X[,i-1] + sqrt(vol[,i-1])*sqrt(dt)*NS
    
    x = vol[,i-1] + kappa*(theta - vol[,i-1])*dt
    y = sqrt(  log(  (xi^2*vol[,i-1]*dt)/(x^2) + 1  )  )
    
    vol[, i] = x*exp(-0.5*y^2+y*NV  )
    
    #Observed Y
    omega = gamma*sqrt(vol[,i])/sqrt(steps)     # n corresponds to steps and not repetitions N? #should vol be i-1? No?
    Y[,i] = X[,i] + omega * rnorm(N, 0, 1)
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

sim.divide <- function(sim.data){
  # divides a list from simulation into Npath sum lists used for pathwise functions (est, Tstat, T* etc)
  # if compact it will not return X or vol
  
  # NEEDS TO BE DONE
  # How do we do this in a smart way?
  #https://stackoverflow.com/questions/12511648/building-a-list-in-a-loop-in-r-getting-item-names-correct ?!
}

