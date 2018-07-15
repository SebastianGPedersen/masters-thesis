setwd(Sys.getenv("masters-thesis"))

### Noise function w. possible AR process
simAR<-function(Npaths, Nsteps, rho, sd){
  
  eps <- matrix(nrow = Npaths, ncol = Nsteps)
  eps[,1]<-rnorm(Npaths, 0, sd)

  for(i in 2:Nsteps){
    eps[,i] <- rho*eps[,i-1]+sqrt(1-rho^2)*rnorm(Npaths,0,sd)
  }
  return(eps)
}

sim.addnoise <- function(process, omega, rho) {
  Npaths <- nrow(process$X)
  Nsteps <- ncol(process$X)
  
  process$Y <- process$X + simAR(Npaths,Nsteps,rho,omega)
  
  return(process)
}

