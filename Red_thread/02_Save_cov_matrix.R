
setwd(Sys.getenv("masters-thesis-data"))
load("BSS_sim/Fit.Rda") #15-20sek. Loads the fit with Gamma function.

setwd(paste0(Sys.getenv("masters-thesis"), "/Vol"))
source("vol_estimators.R")
source("FlexibleFourierForm_Func.R")
source("BV_Analysis_Func.R")
source("PR_Func.R") # PR = Persistence and roughness
source("BSS_model.R") 
source("BSS_Sim.R") 

setwd(paste0(Sys.getenv("masters-thesis"), "/SPY"))
source("datahandling.R")


p0 <- Sys.time() #7sec w. 1000 steps, 1000 paths
Nsteps <- 23400
Npaths <- 1

#Create time points
timepoints <- sim.BSS.equidist_times(Nsteps = Nsteps)

#Save covariance matrix
save.BSS.cov(hVec = timepoints, nPaths = Npaths, S0 = 200, mu_add = 0, type = "Gamma", Fit = Fit)
p0 <- Sys.time() #7sec w. 1000 steps, 1000 paths

