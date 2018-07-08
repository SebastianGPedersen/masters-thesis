setwd(Sys.getenv("masters-thesis"))
source("Module/table1functions.R")
source("Simulation/add_all.R")


####### ESTIMATION PARAMETERS
heston_params <- sim.setup()
h_list <- c(120, 300, 600)/(52*7*24*60*60)
ratio_list <- c(15,12,10) #Same order as h_list
lag <- 10

#Burn a single mu in:
n <- heston_params$Nsteps
mat <- heston_params$mat
dt <- mat/n
n_burn <- max(h_list) / dt
t.index <- seq(from = n_burn, to = 23400, by = 5) #Burn a volatility bandwidth (note 10 in Christensen)

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

c_2_func <- function(beta) {
  c_2 <- 0 #if beta = 0
  if (beta == 0.1) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.25*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.2) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.5*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.3) {
    c_2 <- sqrt((1-2*beta)*(0.00093*0.75*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  else if (beta == 0.4) {
    c_2 <- sqrt((1-2*beta)*(0.00093*1*8)^2/(10/(60*24*7*52))^(1-2*beta))
  }
  return(c_2)
}

alphas <- c(0,0.55,0.65,0.75) #Size and whether or not it is a jump
betas <- c(0,0.1,0.2,0.3,0.4)

#Get burst settings as a list (uses sim.burstsetting for standardization)
burstsettings <- list(length(alphas)*length(betas))
for (beta_index in 1:length(betas)) {
  #beta_index <- 1
  for (alpha_index in 1:length(alphas)){
    #alpha_index <- 1
    burstsettings[[(beta_index-1)*length(alphas)+alpha_index]] <- 
      sim.burstsetting(jump = F, 
                       alpha = alphas[[alpha_index]],
                       beta = betas[[beta_index]], 
                       c_1 = c_1_func(alphas[[alpha_index]]), 
                       c_2 = c_2_func(betas[[beta_index]]))
  }
}

jumps <- c(0.55,0.65,0.75)
for (jump in 1:length(jumps)) {
  burstsettings[[length(burstsettings) + 1]] <- 
    sim.burstsetting(jump = T,
                     alpha = jumps[jump],
                     beta = 0,
                     c_1 = c_1_func(jumps[jump]),
                     c_2 = 0)
}


#Evt. reverse
#### LOOP BECAUSE OF LACK OF MEMORY
Npaths <- 500 #Takes approx. a second per path (because it has to estimate T for 35 processes w. 3 different bandwidths)
n_loops <- ceiling(Npaths/50) #After 50 it just scales linearly if not slower
output_list <- list()
set.seed(100)

p0 <- Sys.time()
for (memory in 1:n_loops) {
    #memory <- 1
  
    ### Keep track
    print(paste("Loop",memory, "out of ",n_loops))
    print(paste("Expected time left:", round(Npaths-(memory-1)*Npaths/n_loops,0),"seconds")) #One second per path
    
    ### SIMULATE ALL PATHS
    heston_params$Npath <- ceiling(Npaths/n_loops)
    Heston <- sim.heston(heston_params)
    
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
Table1 <- restructure_table1(Table1_results,h_list)

### Save and view
save(Table1, file = "Figures2/Saved_data_for_plots/14_table2.Rda")

View(Table1)


