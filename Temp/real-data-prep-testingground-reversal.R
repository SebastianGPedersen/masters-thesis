setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
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
X <- c(600, 300)
Y <- 300
hd <- 1
hv <- 1
t.freq = 5 # every 5 seconds
lag = 10
offset = 10
offset_perId = T
conf = 0.95

id_name = "day"
T_star_TF_col <- 7

# IDENTIFY BURSTS
dt <- data[day <= 3]

# SETUP T-calc every xs per day
grid_T <- data.TforId(data = dt, id_name = id_name, id_number = "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = offset_perId)
grid <- data.TtoStar(grid_T, id_name, conf)

# FUNCTION START FOR REAL
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

test <- reverse_ID(dt, grid_T, grid, 7,
                   X = X, Y = 0, id_name = id_name)

# CHECK THAT IT MAKES SENSE
plot(test$Time)

# CLEAN UP #
# Ideally, we return something like a matrix or clean data table
# Issue is that we want the number of x and ys to be dynamic,
# it is therefore tricky to get table tag as columns from byBurst ( JUST USE PASTE0 ? )

temp<-test[day <= 3,]

# CLEANUP FUNCTION STARTS HERE
dt <- temp

reverse_matrix<-function(reverse_data){
  # CREATES MATRICES FROM REVERSE_DATA
  dt<-reverse_data
  # CREATE X MATRIX
  burst <- dt$burst_no
  nbursts <- length(unique(burst))
  example <- dt[burst == 1,]
  nx <- length(  example[tag == "x", ]$logPrice   )
  ny <- length(  example[ tag == "y", ]$logPrice  )
  
  # SETUP FRAMEWORK
  x <- matrix(NA, nrow = nbursts, ncol = nx)
  xt <-matrix(NA, nrow = nbursts, ncol = nx)
  y <- matrix(NA, nrow = nbursts, ncol = ny)
  #yt <- matrix(NA, nrow = nbursts, ncol = ny)
  ct <- matrix(NA, nrow = nbursts, ncol = 1)
  time <- matrix(NA, nrow = nbursts, ncol = nx+1+ny)
  
  xnames   <- paste0("x", paste0(1:nx))
  ynames   <- paste0("y", paste0(1:ny))
  xtnames  <- paste0("xt", paste0(1:nx))
  #ytnames  <- paste0("yt", paste0(1:ny))
  timename <- c(xnames, "center", ynames)
  
  colnames(x) <- xnames
  colnames(xt) <- xtnames
  colnames(y) <- ynames
  #colnames(yt) <- ytnames
  colnames(ct) <- "ct"
  
  # FILL UP THE MATRICES
  for (i in unique(burst)){
    dat <- dt[burst_no == i,]
    x[i, ] <- dat[tag == "center",]$logPrice - dat[tag == "x",]$logPrice
    y[i, ] <- dat[tag =="y",]$logPrice - dat[tag == "center",]$logPrice
    xt[i, ] <- dat[tag == "x",]$Tval
  #  yt[i, ] <- dat[tag == "y",]$Tval
    ct[i, ] <- dat[tag == "center",]$Tval
    time[i, ] <- dat$DateTime
  }
  return(list(x = x, y = y, xt = xt, ct = ct, time = time))
}

reverse_matrix(test)
# JATAK
