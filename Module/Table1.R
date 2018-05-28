setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")

# INITIAL HESTON SETUP

hestonsetup <- sim.setup(Npath = 10)

# (ONE COULD LOOP OVER MEMORY HERE by setting Npath to something)

# ESTIMATION PARAMETERS
hset <- c(120, 300, 600)/(7*52*24*60*60)
ratio <- 15
lag <- 10
t.index <- seq(from = 2000, to = 23400, by = 5) # tune the initial burn

# SETUP BURST SETTINGS (a = 0, 0.55, 0.65, 0.75 | b = 0.0, 0.1, 0.2, 0.3, 0.4)
burstset11 <- sim.burstsetting(alpha = 0.00, beta = 0.0, c_1 = c_1_func(0.0),  c_2 = c_2_func(0.0))
burstset12 <- sim.burstsetting(alpha = 0.55, beta = 0.0, c_1 = c_1_func(0.55), c_2 = c_2_func(0.0))
burstset13 <- sim.burstsetting(alpha = 0.65, beta = 0.0, c_1 = c_1_func(0.65), c_2 = c_2_func(0.0))
burstset14 <- sim.burstsetting(alpha = 0.75, beta = 0.0, c_1 = c_1_func(0.75), c_2 = c_2_func(0.0))

burstset21 <- sim.burstsetting(alpha = 0.00, beta = 0.1, c_1 = c_1_func(0.0),  c_2 = c_2_func(0.1))
burstset22 <- sim.burstsetting(alpha = 0.55, beta = 0.1, c_1 = c_1_func(0.55), c_2 = c_2_func(0.1))
burstset23 <- sim.burstsetting(alpha = 0.65, beta = 0.1, c_1 = c_1_func(0.65), c_2 = c_2_func(0.1))
burstset24 <- sim.burstsetting(alpha = 0.75, beta = 0.1, c_1 = c_1_func(0.75), c_2 = c_2_func(0.1))

burstset31 <- sim.burstsetting(alpha = 0.00, beta = 0.2, c_1 = c_1_func(0.0),  c_2 = c_2_func(0.2))
burstset32 <- sim.burstsetting(alpha = 0.55, beta = 0.2, c_1 = c_1_func(0.55), c_2 = c_2_func(0.2))
burstset33 <- sim.burstsetting(alpha = 0.65, beta = 0.2, c_1 = c_1_func(0.65), c_2 = c_2_func(0.2))
burstset34 <- sim.burstsetting(alpha = 0.75, beta = 0.2, c_1 = c_1_func(0.75), c_2 = c_2_func(0.2))

burstset41 <- sim.burstsetting(alpha = 0.00, beta = 0.3, c_1 = c_1_func(0.0),  c_2 = c_2_func(0.3))
burstset42 <- sim.burstsetting(alpha = 0.55, beta = 0.3, c_1 = c_1_func(0.55), c_2 = c_2_func(0.3))
burstset43 <- sim.burstsetting(alpha = 0.65, beta = 0.3, c_1 = c_1_func(0.65), c_2 = c_2_func(0.3))
burstset44 <- sim.burstsetting(alpha = 0.75, beta = 0.3, c_1 = c_1_func(0.75), c_2 = c_2_func(0.3))

burstset51 <- sim.burstsetting(alpha = 0.00, beta = 0.4, c_1 = c_1_func(0.0),  c_2 = c_2_func(0.4))
burstset52 <- sim.burstsetting(alpha = 0.55, beta = 0.4, c_1 = c_1_func(0.55), c_2 = c_2_func(0.4))
burstset53 <- sim.burstsetting(alpha = 0.65, beta = 0.4, c_1 = c_1_func(0.65), c_2 = c_2_func(0.4))
burstset54 <- sim.burstsetting(alpha = 0.75, beta = 0.4, c_1 = c_1_func(0.75), c_2 = c_2_func(0.4))


# COMBINE FUNCTIONS AND SETUPS
# This corresponds to all scenarios with the first burstsettings
# (raw, vb, db, vbdb, j, vbj)
# DO THIS ROW-WISE
funs <- list(raw, sim.adddb_arglist, sim.adddb_arglist, sim.adddb_arglist, sim.addjump_arglist, sim.addjump_arglist, sim.addjump_arglist, # row one
             sim.addvb_arglist, vbdb, vbdb, vbdb, vbJ, vbJ, vbJ, # row 2
             sim.addvb_arglist, vbdb, vbdb, vbdb, vbJ, vbJ, vbJ, # row 3
             sim.addvb_arglist, vbdb, vbdb, vbdb, vbJ, vbJ, vbJ, # row 4
             sim.addvb_arglist, vbdb, vbdb, vbdb, vbJ, vbJ, vbJ) # row 5

args <- list(burstset11, burstset12, burstset13, burstset14, burstset12, burstset13, burstset14,  # row 1
             burstset21, burstset22, burstset23, burstset24, burstset22, burstset23, burstset24,  # row 2
             burstset31, burstset32, burstset33, burstset34, burstset32, burstset33, burstset34,  # row 3
             burstset41, burstset42, burstset43, burstset44, burstset42, burstset43, burstset44,  # row 4
             burstset51, burstset52, burstset53, burstset54, burstset52, burstset53, burstset54)  # row 5

# START SIMULATION
sim<-Table1(hest_setup = hestonsetup, fun = funs , args = args,
                       h_list = hset, ratio = ratio, t.index = t.index, lag = lag, conf = 95)
print(sim)
