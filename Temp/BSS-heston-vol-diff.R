setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/BlackScholes.R")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nsteps <- 20000; Npath <- 1000 # 40 000; 10 000 # TAKES ROUGHLY 1.5 hours for all schemes

#hest <- sim.heston(sim.setup(Nsteps = Nsteps, Npath = Npath, mat = 6.5/(6.5*252)))
hest <- sim.heston(sim.setup(Nsteps = Nsteps, Npath = Npath, theta = mean(IV_bss)))

# IMPORT BSS as "HESTON"
setwd(Sys.getenv("masters-thesis-data"))
load("BSSsim.Rdata") # called BSSsim
setwd(Sys.getenv("masters-thesis"))

dt <- (hest$time[2]-hest$time[1])
IV_heston <- rowSums(hest$vol*dt)

dt <- (BSSsim$time[2]-BSSsim$time[1])
IV_bss <- rowSums(BSSsim$vol*dt)

mean(IV_heston)
mean(IV_bss)

sqrt(mean(IV_bss))
