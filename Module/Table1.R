setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")
source("Simulation/add_all.R")


####### ESTIMATION PARAMETERS
heston_params <- sim.setup()
h_list <- c(120, 300, 600)/(52*7*24*60*60)
ratio <- 15
lag <- 100
#Burn a single mu in:
n <- heston_params$Nsteps
mat <- heston_params$mat
dt <- mat/n
n_burn <- max(hset) / dt
t.index <- seq(from = n_burn, to = 23400, by = 5)

# BURST SETTINGS (a = 0, 0.55, 0.65, 0.75 | b = 0.0, 0.1, 0.2, 0.3, 0.4), c_1 and c_2 from Batman
c_1 <- 0.3
c_2 <- 0.016
alphas <- list(c(0,F),c(0.55,T),c(0.55,F),c(0.65,T),c(0.65,F),c(0.75,T),c(0.75,F)) #Size and whether or not it is a jump
betas <- c(0,0.1,0.2,0.3,0.4)

#Get burst settings as a list (uses sim.burstsetting for standardization)
burstsettings <- list(length(alphas)*length(betas))
for (beta_index in 1:length(betas)) {
  for (alpha_index in 1:length(alphas)){
    burstsettings[[(beta_index-1)*length(alphas)+alpha_index]] <- 
      sim.burstsetting(jump = alphas[[alpha_index]][2], alpha = alphas[[alpha_index]][1], beta = betas[[beta_index]][1], c_1 = c_1, c_2 = c_2)
  }
}

#### LOOP BECAUSE OF LACK OF MEMORY
Npaths <- 40 #Takes approx. a second per path (because it has to estimate T for 35 processes w. 3 different bandwidths)
n_loops <- ceiling(Npaths/10) #After 50 it just scales linearly if not slower
output_list <- list()

for (memory in 1:n_loops) {
  
    ### Keep track
    print(paste("Loop",memory, "out of ",n_loops))
    print(paste("Expected time left:", round(Npaths-(memory-1)*Npaths/n_loops,0),"seconds")) #One second per path
    
    ### SIMULATE ALL PATHS
    heston_params$Npath <- ceiling(Npaths/n_loops)
    Heston <- sim.heston(heston_params)
    
    all_simulations <- sim.add_all(Heston = Heston, burst_args = burstsettings)
    
    ### CALCULATE TABLE 1
    output_list[[memory]] <- Table1(all_simulations, h_list = h_list, ratio = ratio, t.index = t.index, lag = lag, conf = 95)
}

### Take mean accross memories (assuming )
Table1_results <- output_list[[1]]
for (memory in 2:n_loops) {
  Table1_results[,3+1:length(h_list)] <- Table1_results[,3+1:length(h_list)] + output_list[[memory-1]][,3+1:length(h_list)]
}
Table1_results[,3+1:length(h_list)] <- Table1_results[,3+1:length(h_list)] /n_loops


###Re-structure Table as in Christensen et. al.:
restructured_table <- restructure_table1(Table1_results,h_list)

### Save and view
#save(resheaped_data, file = "Module/Tableuno.Rdata")
#load("Module/Tableuno.Rdata")
View(restructured_table)


