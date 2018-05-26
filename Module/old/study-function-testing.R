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
h = 300/(3600*24*7*52)
setting<-sim.setup(Nsteps = 23400, Npath = 1000, omega = 1.6*10^-5)
burstset <- sim.burstsetting(alpha = 0.8, beta = 0.1, c_1 = -0.009, c_2 = 0.069)

# find T points
tind<-seq(from = 1000, to = 23399, by = 60)

master <- study(setting = setting, hd = h, hv = h, t.index = tind,
                   conf = 0.990, burstsetting = burstset)

master$output

if(F){
  # DEBUG TEST #
  data<-master$debug
  data$Y <- t(diff(t(data$Y)))
  
  mu<-est.mu.mat.next(data, h, tind)
  sig <- est.sigma.mat.next(data, h, tind, lag = 10)
  
  test<-sqrt(h/0.5)*mu$mu/sqrt(sig$sig)
  
  Tstar <- apply(abs(test), 1, max)
  master$output$Tstarraw-Tstar
  
  plot(master$debug$Y[1,], type="l")
  
  # DEBUG SIM #
  alpha <- burstset$alpha;  beta <- burstset$beta;
  burst_time <- burstset$burst_time;  interval_length <- burstset$interval_length;
  c_1 <- burstset$c_1;  c_2 <- burstset$c_2
}

# abe <- function(data) addDB(data = data, amroariogjarg)



