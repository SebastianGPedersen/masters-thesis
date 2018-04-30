#######################################################################
#                                                                     #
# FROM HERE WE CHAIN EVERYTHING TOGETHER TO FORM THE SIMULATION STUDY #
#                                                                     #
#######################################################################

# working directory should be "masters-thesis"
setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
source("module/SimStudyFunction.R")

# params
hset = c(120,300,600)/(3600*24*7*52)
alphaset = c(0.55, 0.65, 0.75)
betaset = c(0, 0.1, 0.2 ,0.3 ,0.4)
setting<-sim.setup(Nsteps = 23400, Npath = 1)

# find T points
tind<-seq(from = 60, to = 23399, by = 60)

masterLength <- length(betaset)*length(alphaset)*length(hset)
master <- numeric(masterLength) 
i <- 0
timeElapsed <- proc.time()[3]

for(beta in betaset){
  for(alpha in alphaset){
    for(h in hset){
      i <- i + 1
      # if(i>10){
      timeElapsed <-  (proc.time()[3]- timeElapsed)
        print(c(i, masterLength, timeElapsed))
        
        #create setting
        burstset <- sim.burstsetting(alpha = alpha, beta = beta,
                                     burst_time = burst_time, interval_length = interval_length,
                                     c_1 = 0.1, c_2 = 0.1)
        
        # master will be a vector of lists(h/alpha/beta) of lists(input/output)
        master[i] <- study(setting = setting, hd = h, hv = h, t.index = tind,
                        conf = 0.95, burstsetting = burstset)
      # }
      
    }
  }
}


