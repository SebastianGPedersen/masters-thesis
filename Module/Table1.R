setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")

# INITIAL HESTON SETUP

hestonsetup <- sim.setup(Npath = 10)

# (ONE COULD LOOP OVER MEMORY HERE by setting Npath to something)

# ESTIMATION PARAMETERS
hset <- c(120, 300, 600)/(7*52*24*60*60)
ratio <- 15
lag <- 10
t.index <- seq(from = 1000, to = 23400, by = 5) # tune the initial burn

# SETUP BURST SETTINGS
burstset11 <- sim.burstsetting(alpha = 0.55, beta = 0.1, c_1 = c_1_func(0.55), c_2 = c_2_func(0.1))
burstset12 <- sim.burstsetting(alpha = 0.55, beta = 0.2, c_1 = c_1_func(0.55), c_2 = c_2_func(0.2))
burstset13 <- sim.burstsetting(alpha = 0.55, beta = 0.3, c_1 = c_1_func(0.55), c_2 = c_2_func(0.3))
burstset14 <- sim.burstsetting(alpha = 0.55, beta = 0.4, c_1 = c_1_func(0.55), c_2 = c_2_func(0.4))

# COMBINE FUNCTIONS AND SETUPS
# This corresponds to all scenarios with the first burstsettings
# (raw, vb, db, vbdb, j, vbj)
funs <- list(raw, sim.addvb_arglist, sim.adddb_arglist, vbdb, sim.addjump_arglist, vbJ)
args <- list(burstset11, burstset11, burstset11, burstset11, burstset11, burstset11)

# START SIMULATION
sim<-Table1(hest_setup = hestonsetup, fun = funs , args = args,
                       h_list = hset, ratio = ratio, t.index = t.index, lag = lag, conf = 95)
print(sim)
