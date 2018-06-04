reverse_ID <- function(data, T_grid_data, T_star_data, T_star_TF_col,
                       X, Y, id_name){
  # INPUT (data) SHOULD BE MARKED IN DAYS
  dt <- data
  grid_T <- T_grid_data
  grid <- T_star_data
  
  # dates for matching
  dates <- dt[["DateTime"]]
  dates_T <- grid_T[["DateTime"]]
  
  perID <- function(by, index){
    # ID DAY STARTS HERE
    id <- by
    # DEFINE BOUNDARY for each id
    l_bound <- id[1, DateTime]; u_bound <- id[length(id[,DateTime]), DateTime];
    burst <- grid[grid[[T_star_TF_col]]== T, ]
    
    perBurst<-function(index){
      # LETS CHOOSE A POINT (WE WOULD DO THIS BY ID AGAIN)
      center <- burst[burst[[id_name]] == index,]
      
      prev <- center$DateTime-X
      after <- center$DateTime+Y
      
      # GET DATA
      if(any(prev < l_bound || after > u_bound)){
        prev_index <- rep(NA, length(prev))
        after_index <- rep(NA, length(after))
        return() # we wish to return empty table or NAs
      }
      else{
        # match dates roughly
        prev_index <- data.date_To_tindex(dates = prev, dates)
        after_index <- data.date_To_tindex(dates = after, dates)
      }
      
      # PACK UP DATA NICELY
      prev_data<-dt[prev_index, ]
      after_data <- dt[after_index, ]
      
      # setup center data
      center_index <- match(center$DateTime, dates)
      center_data <- dt[center_index, ]
      
      # CALC T
      # FIND PREVMU/SIG
      Tprev_index <- sum(grid_T[["DateTime"]] < min(prev))
      
      prevmu <- list(time = grid_T[Tprev_index,]$Time, mu = grid_T[Tprev_index,]$Mu)
      prevsig <- list(time = grid_T[Tprev_index,]$Time, sig = grid_T[Tprev_index,]$Sigma)
      
      # find sub data for estimation
      estsub <- dt[Tprev_index:max(prev_index),]
      estdata_index <- match(prev_data$Time,estsub$Time)
      estdata <- list(time = estsub$Time, Y = diff(estsub$logPrice))
      
      #estimation
      mu <- est.mu.next(estdata, prevmu, hd, estdata_index)
      sig <- est.sigma.next(estdata, prevsig, hv, estdata_index, lag = lag) # DEN SKAL IKKE RETURNE PREVSIG
      testT <- teststat(mu, sig, hd, hv)$test
      
      prev_data[, Tval := testT]
      center_data[, Tval := center$Tstar]
      after_data[, Tval := 0] # We could calc this value
      
      # TAG 'EM
      prev_data <- prev_data[, tag := "x"]
      after_data <- after_data[, tag:= "y"]
      center_data <- center_data[, tag := "center"]
      # PAK VARERNE
      out <- rbind(prev_data, after_data, center_data )
      return(out)
    }
    
    # only those with bursts
    out<-burst[, perBurst(.BY), by = burst[[id_name]]]
  }
  out<-dt[, perID(.SD, .BY), by = dt[[id_name]]]
  # MARK THEM AFTER IT HAS BEEN DONE FOR EACH DAY
  nbursts <- length(out[tag == "center"]$Time)
  entry_length <- length(out$Time)/nbursts
  burst_no <- rep(1:length(out[tag == "center"]$Time), each = entry_length)
  out <- out[, burst_no := burst_no]
  out$burst <- NULL
  return(out)
}