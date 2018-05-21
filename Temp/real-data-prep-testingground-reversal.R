setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
source("module/SimStudyFunction.R")
source("spy/dataFunctions.R")
source("spy/datahandling.R")

# LET US TAKE A LOOK AT SECTION 5.3 - REVERSION AFTER DRIFT BURSTS

# The idea is that we can 'model' the return up to X-min after the burst peak
# by looking at the return X-min down to Y-min before
# Ideally, we could expand this into something more elaborate and ~cool~
# We CANNOT use T in the model unless we calculate T way more often... (Once per second)
# OR we could calculate Ts whenever we need them - on the fly inside function...?

# FETCH DATA AND LABEL DAYS
# by labeling days, we 
data<-data.getFull()
data<-data.dayID(data, id = "day")


#FUNCTION START

sub_id <- "5min"
X <- c(600, 300)
Y <- 300
gridwidth <- 300
hd <- 1
hv <- 1
t.freq = 5
lag = 10
offset = 10
offset_perId = F
conf = 0.95

# We would use this function by id - we extract a single day to simulate the by id effect.
if(F){
  # ID DAY STARTS HERE
  dt <- data[day == 3,]
  
  # DEFINE BOUNDARY (9:30 - 15:30) for each day
  l_bound <- dt[1, DateTime]; u_bound <- dt[length(dt[,DateTime]), DateTime];
  
  # THE GRID - A DIGITAL FRONTIER - I TRIED TO PICTURE CLUSTERS OF INFORMATION AS THEY MOVED THROUGH THE COMPUTERS
  # We divide every day into several X-sec intervals
  subid <- data.xsecID(dt, gridwidth, sub_id)
  grid_T <- data.TforId(subid, sub_id, hd = hd, t.freq = t.freq, lag = lag, offset = offset, offset_perId = offset_perId) #UPDATE WITH HV
  grid <- data.TtoStar(grid_T, sub_id, conf = conf)
  
  dates <- subid[["DateTime"]]
  dates_T <- grid_T[["DateTime"]]
  # Choose only cases with burst (col 7 is the T/F column)
  burst <- grid[grid[[7]]== T, ] # This only works on standard data due to hardcoded index
  
  # We mark each of the bursts with an ID 'burst_no'
  burst <-burst[, burst_no := 1:length(burst$DateTime)]
  
  # ID 'SUBID' STARTS HERE
  perBurst<-function(index){
    # LETS CHOOSE A POINT (WE WOULD DO THIS BY ID AGAIN)
    center <- burst[burst[[sub_id]] == index,]
    burst_index <- index
    
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
    prev_data<-subid[prev_index, ]
    after_data <- subid[after_index, ]
    
    # setup center data
    center_index <- match(center$DateTime, dates)
    center_data <- subid[center_index, ]
    
    # CALC T
    # FIND PREVMU/SIG
    Tprev_index <- sum(grid_T[["DateTime"]] < min(prev))
  
    prevmu <- list(time = grid_T[Tprev_index,]$Time, mu = grid_T[Tprev_index,]$Mu)
    prevsig <- list(time = grid_T[Tprev_index,]$Time, sig = grid_T[Tprev_index,]$Sigma)
    
    estsub <- subid[Tprev_index:max(prev_index),]
    estdata_index <- match(prev_data$Time,estsub$Time)
    estdata <- list(time = estsub$Time, Y = diff(estsub$logPrice))
    
    mu <- est.mu.next(estdata, prevmu, hd, estdata_index)
    sig <- est.sigma.next(estdata, prevsig, hv, estdata_index, lag = lag) # DEN SKAL IKKE RETURNE PREVSIG
    
    testT <- teststat(mu, sig, hd, hd)$test
    
    prev_data[, Tval := testT]
    center_data[, Tval := center$Tstar]
    after_data[, Tval := 0] # We could calc this value
    
    
    # TAG 'EM
    prev_data <- prev_data[, tag := "x"]
    after_data <- after_data[, tag:= "y"]
    center_data <- center_data[, tag := "center"]
    out <- rbind(prev_data, after_data, center_data)
    
    #Mark by index
    #if(any(names(out)=="id")) stop("Data table id already exists!")
    #out[, id:=index]
    
    return(out)
  }
  
  sub<-burst[, perBurst(.BY), by = burst[[sub_id]]]
  
}

perDay <- function(data, index){
  # ID DAY STARTS HERE
  dt <- data
  # DEFINE BOUNDARY for each day
  l_bound <- dt[1, DateTime]; u_bound <- dt[length(dt[,DateTime]), DateTime];
  
  # THE GRID - A DIGITAL FRONTIER - I TRIED TO PICTURE CLUSTERS OF INFORMATION AS THEY MOVED THROUGH THE COMPUTERS
  # We divide every day into several X-sec intervals
  subid <- data.xsecID(dt, gridwidth, sub_id)
  grid_T <- data.TforId(subid, sub_id, hd = hd, t.freq = t.freq, lag = lag, offset = offset, offset_perId = offset_perId) #UPDATE WITH HV
  grid <- data.TtoStar(grid_T, sub_id, conf = conf)
  
  dates <- subid[["DateTime"]]
  dates_T <- grid_T[["DateTime"]]
  # Choose only cases with burst (col 7 is the T/F column)
  burst <- grid[grid[[7]]== T, ] # This only works on standard data due to hardcoded index
  
  # We mark each of the bursts with an ID 'burst_no'
  burst <-burst[, burst_no := 1:length(burst$DateTime)]
  
  # ID 'SUBID' STARTS HERE
  perBurst<-function(index){
    # LETS CHOOSE A POINT (WE WOULD DO THIS BY ID AGAIN)
    center <- burst[burst[[sub_id]] == index,]
    burst_index <- index
    
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
    prev_data<-subid[prev_index, ]
    after_data <- subid[after_index, ]
    
    # setup center data
    center_index <- match(center$DateTime, dates)
    center_data <- subid[center_index, ]
    
    # CALC T
    # FIND PREVMU/SIG
    Tprev_index <- sum(grid_T[["DateTime"]] < min(prev))
    
    prevmu <- list(time = grid_T[Tprev_index,]$Time, mu = grid_T[Tprev_index,]$Mu)
    prevsig <- list(time = grid_T[Tprev_index,]$Time, sig = grid_T[Tprev_index,]$Sigma)
    
    estsub <- subid[Tprev_index:max(prev_index),]
    estdata_index <- match(prev_data$Time,estsub$Time)
    estdata <- list(time = estsub$Time, Y = diff(estsub$logPrice))
    
    mu <- est.mu.next(estdata, prevmu, hd, estdata_index)
    sig <- est.sigma.next(estdata, prevsig, hv, estdata_index, lag = lag) # DEN SKAL IKKE RETURNE PREVSIG
    
    testT <- teststat(mu, sig, hd, hd)$test
    
    prev_data[, Tval := testT]
    center_data[, Tval := center$Tstar]
    after_data[, Tval := 0] # We could calc this value
    
    
    # TAG 'EM
    prev_data <- prev_data[, tag := "x"]
    after_data <- after_data[, tag:= "y"]
    center_data <- center_data[, tag := "center"]
    out <- rbind(prev_data, after_data, center_data)
    
    #Mark by index
    #if(any(names(out)=="id")) stop("Data table id already exists!")
    #out[, id:=index]
    
    return(out)
  }
  
  sub<-burst[, perBurst(.BY), by = burst[[sub_id]]]
}

dt <- data[data[["day"]] <= 10]

# THIS TAKES TIME
dt[, perDay(.SD, .BY), by = dt[["day"]]]

# CLEAN UP #
# Ideally, we return something like a matrix or clean data table
# Issue is that we want the number of x and ys to be dynamic,
# it is therefore tricky to get table tag as columns from byBurst

# CLEANUP WILL BE TOMORROW
# For each burst, take out all x, y and center and flatten info from each burst into a single row

