#### source functions neeeded

setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")

cd_UNIQUE_NAME_BSS_Sim <- getwd()

setwd(paste0(Sys.getenv("masters-thesis"), "/Vol"))

source("vol_estimators.R")
source("FlexibleFourierForm_Func.R")
source("BV_Analysis_Func.R")
source("PR_Func.R") # PR = Persistence and roughness
source("BSS_model.R") 
source("BSS_Sim.R") 

setwd(paste0(Sys.getenv("masters-thesis"), "/SPY"))
source("datahandling.R")

setwd(cd_UNIQUE_NAME_BSS_Sim)
rm(cd_UNIQUE_NAME_BSS_Sim)
##################################
# temp <- BV.get_SPY_data(maxDate = "2014-01-03")
# plot(temp$DateTime, exp(temp$logPrice))

Fit <- sim.BSS.Fit(bucketLengthInMinutes = 5)
Fit$bvS_List$bvSDTfff
Fit$alpha
Fit$memory_param
Fit$log_c_sigma
Fit$nu

save(Fit, file = "Fit.Rdata")

setting <- sim.setup(Npath = 2, Nsteps = 23400)
timepoints <- sim.heston.uneven(setting)$time

uneven <- sim.BSS(hVec = timepoints, nPaths = 1000, S0 = 200, mu_add = 0, type = "Gamma", Fit = Fit)

### Save and view
save(uneven, file = "Uneven.Rdata")
