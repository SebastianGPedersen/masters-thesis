setwd(Sys.getenv("masters-thesis"))
source("Module/Table1functions.R")
source("Simulation/add_all.R")


####### ESTIMATION PARAMETERS
h_list <- c(120, 300, 600)/(52*7*24*60*60)
ratio_list <- c(15,12,10)
lag <- 10

# IMPORT BSS as "HESTON"
setwd(Sys.getenv("masters-thesis-data"))
load("BSSsim.Rdata") # called BSSsim
setwd(Sys.getenv("masters-thesis"))

setting <- sim.setup(Nstep = 20000)

BSSsim$Y <- BSSsim$X + setting$omega*rnorm(   length(BSSsim$time), 0 , 1   ) # add microstruct

# Burn a single mu in:
n <- length(BSSsim$time)
mat <- 6.5/(52*7*24)
dt <- mat/n
n_burn <- ceiling(max(h_list) / dt)
t.index <- seq(from = n_burn, to = n, by = 5) #Burn a volatility bandwidth (note 10 in Christensen)

# BURST SETTINGS (a = 0, 0.55, 0.65, 0.75 | b = 0.0, 0.1, 0.2, 0.3, 0.4), c_1 and c_2 from Batman
c_1_func <- function(alpha) {
  c_1 <- 0 #If alpha == 0
  if (alpha == 0.55) {
    c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)
  } 
  else if (alpha == 0.65){
    c_1 <- (1-alpha)*0.01/(10/(60*24*7*52))^(1-alpha)
  } 
  else if (alpha == 0.75) {
    c_1 <- (1-alpha)*0.015/(10/(60*24*7*52))^(1-alpha)
  }
  return(c_1)
}

c_2_func <- function(beta) { # multiply 8 infront of 0.00093 for high vol #
  c_2 <- 0 #if beta = 0
  if (beta == 0.1) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.25)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.2) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.5)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.3) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.75)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.4) {
    c_2 <- sqrt((1-2*beta)*(0.00093*1)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  return(c_2)
}

alphas <- list(c(0,F),c(0.55,F),c(0.55,T),c(0.65,F),c(0.75,F)) #Size and whether or not it is a jump
betas <- c(0,0.1,0.2,0.3,0.4)

#Get burst settings as a list (uses sim.burstsetting for standardization)
burstsettings <- list(length(alphas)*length(betas))
for (beta_index in 1:length(betas)) {
  #beta_index <- 1
  for (alpha_index in 1:length(alphas)){
    #alpha_index <- 1
    burstsettings[[(beta_index-1)*length(alphas)+alpha_index]] <- 
      sim.burstsetting(jump = alphas[[alpha_index]][2], 
                       alpha = alphas[[alpha_index]][1], 
                       #reverse = T,
                       #recenter = T,
                       beta = betas[[beta_index]][1], 
                       c_1 = c_1_func(alphas[[alpha_index]][1]), 
                       c_2 = c_2_func(betas[[beta_index]][1]))
  }
}

#Evt. reverse
#### LOOP BECAUSE OF LACK OF MEMORY
Npaths <- 1000 #Takes approx. a second per path (because it has to estimate T for 35 processes w. 3 different bandwidths)
n_loops <- ceiling(Npaths/50) #After 50 it just scales linearly if not slower
output_list <- list()

pathsPer <- Npaths/n_loops

p0 <- Sys.time()
for (memory in 1:n_loops) {
  #memory <- 1
  
  ### Keep track
  print(paste("Loop",memory, "out of ",n_loops))
  print(paste("Expected time left:", round(Npaths-(memory-1)*Npaths/n_loops,0),"seconds")) #One second per path
  
  ### INSTEAD OF SIMULATION WE EXTRACT WHAT WE NEED
  Heston <- list(time = BSSsim$time,
                 Y = BSSsim$Y[((memory-1)*pathsPer+1):(memory*pathsPer), ],
                 X = BSSsim$X[((memory-1)*pathsPer+1):(memory*pathsPer), ],
                 vol = BSSsim$vol[((memory-1)*pathsPer+1):(memory*pathsPer), ]
                 )
  
  all_simulations <- sim.add_all(Heston = Heston, burst_args = burstsettings)
  
  ### CALCULATE TABLE 1
  output_list[[memory]] <- Table1_func(all_simulations, h_list = h_list, ratio = ratio_list, t.index = t.index, lag = lag, conf = 95)
}
print(Sys.time()-p0)

### Take mean accross memories (assuming same number of paths in every memory loop)
Table1_results <- output_list[[1]]
for (memory in 2:n_loops) {
  Table1_results[,3+1:length(h_list)] <- Table1_results[,3+1:length(h_list)] + output_list[[memory-1]][,3+1:length(h_list)]
}
Table1_results[,3+1:length(h_list)] <- Table1_results[,3+1:length(h_list)] /n_loops


###Re-structure Table as in Christensen et. al.:
Table2 <- restructure_table1(Table1_results,h_list)

### Save and view
save(Table2, file = "Module/BSS_Table2_normalvol.Rdata")

View(Table2)


