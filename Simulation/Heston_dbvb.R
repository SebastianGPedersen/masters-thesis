source(paste(Sys.getenv("masters-thesis"),"Simulation/Heston.R",sep="/"))


sim.all <- function(settings,t_burst,t_interval,alpha,beta,c_1,c_2) {
  
  #Initializations
  time = 1:settings$steps
  
  X = matrix(nrow = settings$Npaths, ncol = settings$steps)
  Y = matrix(nrow = settings$Npaths, ncol = settings$steps)
  vol = matrix(nrow = settings$N, ncol = settings$steps)
  
  #Define burst interval
  burst_begin = t_burst-interval/2
  burst_end = t_burst+interval/2
  mat = settings$mat #save maturity
  steps = settings$Nsteps #save nsteps
  
  if (burst_begin < 0 | burst_end > ttm) {
    stop("Burst interval is not contained in simulation interval")
  }
  
  
  # ----------------------- Heston simulation before burst interval ----------------------
  
    #Arguments for sim.heston (take steps until inside interval)
    settings$Nsteps = ceil((burst_begin/mat)*steps)
    settings$mat = burst_begin
    
    #Points before interval
    (time_b,Y_b,X_b,vol_b) = sim.heston(settings)
  
  
  # ----------------------- Heston, VB and DB+VB inside interval ----------------------
  
    #Arguments for sim.heston
    settings$Nsteps = floor(burst_end/mat*steps+1) - settings$Nsteps #Because of "]", we take another step inside on boundary
    settings$mat = ... - settings$mat
    
    #Points inside interval
    (time_H,Y_H,X_H,vol_H) = sim.heston
    (time_vb,Y_vb,X_vb,vol_vb) = sim.vb
    (time_dbvb,Y_dbvb,X_dbvb,vol_dbvb) = sim.dbvb

  
  # ----------------------- Heston again outside interval ----------------------
    
    #Arguments for sim.heston
    
    
    #Points after interval
    (time_H_a,Y_H_a,X_H_a,vol_H_a) = sim.heston
    (time_vb_a,Y_vb_a,X_vb_a,vol_vb_a) = sim.heston
    (time_dbvb_a,Y_dbvb_a,X_dbvb_a,vol_dbvb_a) = sim.heston 

        
  # ----------------------- Wrap together and return ----------------------
    
}

