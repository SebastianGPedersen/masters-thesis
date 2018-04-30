# Simulation of simple BS

sim.BlackScholes <- function(mean, sd, omega, mat = 6.5/(52*7*24), Nsteps = 23400, Npath = 1000){
  # Simulates BS paths - simple...
  Y <- matrix(NA, Npath, Nsteps)
  dt <- (mat/Nsteps)
  for(i in 1:Npath)
  {
    x <- rnorm(n = Nsteps , mean = mean, sd = sd*sqrt(dt))
    eps <- rnorm(n = Nsteps, mean = 0, sd = omega)
    x <- cumsum(x)
    Y[i,] <- x + eps
  }
  
  return(list(time = 1:Nsteps*dt, Y = Y))
}