setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
source("estimation/estimates_reloaded.R")
source("estimation/estimates_revolution.R")
source("vol/vol_estimators.R")
source("vol/FlexibleFourierForm_Func.R")
source("vol/BV_Analysis_Func.R")
source("vol/PR_Func.R") # PR = Persistence and roughness
source("vol/BSS_model.R") 
source("vol/BSS_Sim.R") 

# CHECK FOR SEASONALITY OF T IN UN-EVEN DT STUDY TO DETERMINE ORIGIN OF SEASONALITY

# SIMULATE
setting <- sim.setup(Npath = 1000, Nsteps = 20000, omega = 1.6*10^-5)

seed<-2342
set.seed(seed)

# IMPORT BSS as "HESTON"
setwd(Sys.getenv("masters-thesis-data"))
load("BSSsim.Rdata") # called BSSsim
setwd(Sys.getenv("masters-thesis"))
load("Fit.Rdata")

BSSsim$Y <- BSSsim$X + setting$omega*rnorm(   length(BSSsim$time), 0 , 1   ) # add microstruct

# LOG RETURNS
BSSsim$Y <- t(diff(t(as.matrix(BSSsim$Y))))

# SPLIT UP DATA
BSS.season <- function(hVec, nPaths, S0 = 1, mu_add = 0, type = "Gamma", Fit){
  time_points <- hVec[-1] #For output
  # Assums hVec[1] = 0
  if(missing(Fit)){
    Fit <- sim.BSS.Fit()
  }
  
  hVec <- hVec * (60*24*7*52)/Fit$bvS_List$bucketLengthInMinutes # time in buckets
  dt   <- diff(hVec)
  hVec <- hVec[-1] #Remove time 0
  
  seasonal_Component <- sim.BSS.Seasonality(hVec, Fit$bvS_List)
  return(list(time = time_points, vol = exp(seasonal_Component)))
}

# ESTIMATION PARAMETERS
hd <- 300/(52*7*24*60*60) #(seconds)
hv <- 12*hd                                          
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)
t.index <- seq(offset, setting$Nsteps, by = t.freq)

# ESTIMATION (OF HALF OF THE DATA)
mu <- sqrt(hd)*est.mu.mat.2.0(data = BSSsim, hd = hd, bandwidth_rescale = T)$mu[1:500, t.index]
sigma2 <- est.sigma.mat.2.0(data = BSSsim, hv = hv, lag = lag, bandwidth_rescale = T)$sig[1:500, t.index]
Tstat <- mu/sqrt(sigma2)

# SEASON
#season <- BSS.season(hVec = BSSsim$time, nPaths = 1000)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MEAN ACROSS DAYS
mus<-colMeans(mu)
sis<-colMeans(sigma2)

data<-data.frame(time = BSSsim$time[t.index], Tstat = meanT, Season = season$vol[t.index])

require(ggplot2)
ggplot() +
  geom_point(data=data, aes(x=time*60*60, y=Tstat, col = "Non-equidistant BSS"), size = 1) +
  geom_point(data=data, aes(x=time*60*60, y=Season, col = "Non-equidistant Heston"), size = 1) +
  xlab("Time") + ylab("Average Sigma estimator") 
