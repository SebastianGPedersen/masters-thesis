require(data.table)

data.getFull <- function(){
  options(digits = 14, digits.secs = 5)
  
  #### source and get data - slow
  cd <- getwd()
  setwd(paste0(Sys.getenv("masters-thesis"),"/SPY"))
  source("dataFunctions.R")
  
  setwd(paste0(Sys.getenv("masters-thesis-data"),"/SPY"))
  fullData<-readRDS("2014_SPY_Vol_Avg.rds")
  setwd(cd)
  rm(cd)
  return(fullData)
}

data.getbitdata <- function(name){
  options(digits = 14, digits.secs = 5)
  
  cd <- getwd()
  #### source and get data - slow
  setwd(paste0(Sys.getenv("masters-thesis"),"/SPY"))
  source("dataFunctions.R")
  
  setwd(paste0(Sys.getenv("masters-thesis-data"),"/Bitcoin"))
  fullData<-readRDS(paste0(name,".rds"))
  setwd(cd)
  rm(cd)
  return(fullData)
}

data.dayID <- function(datatable, id = "day"){
  # datatable should be a data.table.
  # id denotes the name of the id column
  dt <- datatable
  # sort by dateTime
  setkey(dt, DateTime)
  if("Date" %in% names(dt)){
    # Identify unique days
    unDates <- unique(dt$Date)
    # Match with ID number
    dtB<-match(dt$Date, unDates)
  }
  else{
    unDates <- unique(as.Date(dt$DateTime))
    dtB<-match(as.Date(dt$DateTime), unDates)
  }
  
  # Add ID column
  dt[, paste0(id) := dtB]
  gc()
  return(dt)
}

data.xsecID <- function(datatable, bucketLengthInSeconds, id = "id"){
  # datatable should be a data.table
  # id denotes the name of the id column
  dt<-copy(datatable)       # NOT SURE WHY, BUT COPY FIXED A BUG ABOUT .SD BEING LOCKED
  time <- dt[, DateTime]
  #Intraday
  buckets <- seq(data.floor_date(time[1]), time[length(time)], by = bucketLengthInSeconds)
  buckets[1] <- buckets[1] - 0.0000001
  dtB2<- .bincode(dt$DateTime, breaks = c(buckets, time[length(time)]))
  dt[, paste0(id) := dtB2]
  return(dt)
}

data.TforId <- function(data, id_name, id_number, hd, hv, t.freq, lag, offset, offset_perId = T){
  # Estimates T test for every t.freq (unit: seconds) within a certain id
  # hd and lag are parameters for estimation
  # offset as vector, where each element is how many observations are skipped for that id
  # If offset_perId is false, then offset is only applied on first
  # If offset_perId is true, but offset is only a single value, then this value will be used for every id
  # Bandwidth, freq, lag and offset can all be vectors
  
  # IDENTIFY COLUMN INDEX
  if(id_name %in% names(data)){
    id_index <- match(id_name, names(data))
  }
  else{
    stop("id_name not found in data")
  }
  
  if(missing(id_number)){
    id_number <- unique(data[[id_index]])
  }
  if(id_number == "all"){
    id_number <- unique(data[[id_index]])
  }
  # PREVENT USER FROM FKNG UP THE DUMB LOGIC OF OFFSETS
  if(length(offset) > 1 & offset_perId == F){
    stop("Several offsets given (vector) but offset_perId is False - This is a contradiction")
  }
  if(length(offset) == 1){
    if(offset_perId){
      offset <- rep(offset, length(id_number))
    }
    else{
      offset <- c(offset, rep(0, length(id_number)-1))
    }
  }
  if(length(offset) != length(id_number) & offset_perId == TRUE){
    stop("With offset_perId, the length of offset should either be one or equal to length of id_number")
  }
  
  if(length(hd) == 1) hd <- rep(hd, length(id_number))
  if(length(hd) != length(id_number)) stop("Length of hd should either be one or matching length of id_number")
  
  if(length(t.freq) == 1) t.freq <- rep(t.freq, length(id_number))
  if(length(t.freq) != length(id_number)) stop("Length of t.freq should either be one or matching length of id_number")
  
  if(length(lag) == 1) lag <- rep(lag, length(id_number))
  if(length(lag) != length(id_number)) stop("Length of lag should either be one or matching length of id_number")
  
  # ALL SINGLE/VECTOR TROUBLES SHOULD BE OVER BY NOW
  
  # DATA
  dt <- data[data[[id_index]] %in% id_number, ] # this should extract correct columns
  
  ids <- dt[[id_index]]
  # DEFINE INTERNAL FUNCTION FOR HANDLING .BY ID
  TforId.internal<-function(data, index){
    
    if(is.list(index)) index <- index[[1]]
    
    id_index<-match(index, unique(ids))
    # INDEX HANDLING
    hd <- hd[id_index]
    t.freq <- t.freq[id_index]
    lag <- lag[id_index]
    offset <-offset[id_index]
    #
    
    dt<-data
    # data input are relevant IDs
    dy <- dt[, lapply(.SD, diff), .SDcols = "logPrice"]$logPrice
    times <- dt[, .SD, .SDcols = c("DateTime","Time")]
    time <- times$Time; Date <- times$DateTime
    
    # DETERMINE INPUTS
    input <- list(time = time, Y = dy)
    
    # SETUP T-INDEX #
    dates<-seq(data.floor_date(Date[1]), Date[length(Date)], by = t.freq)
    t.index<-data.date_To_tindex(dates, t_dates = Date)
    
    # ENSURE ONLY UNIQUE
    t.index <- unique(t.index)
    
    # OFFSET
    t.index <- t.index[(1+offset):length(t.index)]
    
    # <ESTIMATION>
    mu<-est.mu.next.cpp(input, hd = hd, t.index = t.index)
    sig<-est.sigma.next.cpp(input, hv = hv, t.index = t.index, lag = lag)
    
    # RESCALING
    mu_re <- est.rescale.mu.vec(mu$mu, mu$time, t_beginning = time[1], h_mu = hd)
    sig_re <- est.rescale.sigma.vec(sig$sig, sig$time, t_beginning = time[1], h_sigma = hv)
    
    Ttest <- sqrt(hd)*mu_re/sqrt(sig_re)
    mu <- mu_re
    sig <- sig_re 
    # </ESTIMATION>
    
    # FORMULATE OUTPUT
    Nout<-length(t.index)
    out <- data.table(DateTime = Date[t.index], Time = time[t.index], Mu = mu, Sigma = sig, Tval = Ttest)
    out <- out[, paste0(id_name) := index]
    return(out)
  }
  
  # APPLY FUNCTION
  out<-dt[, TforId.internal(.SD, .BY), by = dt[[id_index]]]
  
  return(out)
}

data.TstarforId <- function(data, id_name, id_number, hd, hv, t.freq, lag, offset, offset_perId = T, conf){
  # WRAPPER AROUND data.TforId THAT RETURNS THE TSTAR VALUE INSTEAD
  if(missing(id_number)){
    Ttest <- data.TforId(data = data, id_name = id_name, id_number = id_number, hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = offset_perId)
  }
  else{
    Ttest <- data.TforId(data = data, id_name = id_name, hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = offset_perId)
  }
  
  Tstar.internal <- function(data, index, conf){
    # UNPACK
    if(is.list(index)) index <- index[[1]]
    
    info <- data[, .SD, .SDcols = c("DateTime", "Tval","Time")]
    dates <- info$DateTime; Ts <- info$Tval; time <- info$Time
    #GET OUT
    Tstar <- as.numeric(max(abs(Ts)))
    t.point <- match(Tstar, abs(Ts))
    # FIND ID
    out <- data.table(DateTime = dates[t.point], Time = time[t.point], Tstar = Tstar)
    out[, paste0(id_name):= index]
    
    if(missing(conf)){
      return(out)
    }
    # IF CONF IS ABOVE ONE - THEN WE ASSUME IT IS A THRESHOLD INSTEAD
    mode <- any(conf >= 1)
    if(!mode){
      # FIT RHO
      rho <- est.rho.MLE(Ts)
      rhom <- length(Ts)
      rhorho <- rho
      
      # THE POSSIBILITY OF MULTIPLE CONFIDENCES (can become messy real quick)
      if(length(conf) > 1){
        for(i in 1:length(conf)){
          quantile <- est.z_quantile(rhom, rhorho, conf[i])$qZm
          out[, paste0("q_",conf[i]):= quantile]
          out[, paste0("db_",conf[i]) := Tstar >= quantile]
          
        }
      }
      else{
        quantile <- est.z_quantile(rhom, rhorho, conf)$qZm
        out[, paste0("q_",conf):= quantile]
        out[, paste0("db_",conf) := Tstar >= quantile]
      }
    }
    else{
      # THE POSSIBILITY OF MULTIPLE CONFIDENCES (can become messy real quick)
      if(length(conf) > 1){
        for(i in 1:length(conf)){
          out[, paste0("threshold_",conf[i]):= conf[i]]
          out[, paste0("db_",conf[i]) := Tstar >= conf[i]]
          
        }
      }
      else{
        out[, paste0("threshold__",conf):= conf]
        out[, paste0("db_",conf) := Tstar >= conf]
      }
    }
    return(out)
    
  }
  # IDENTIFY ID COLUMN IN Ttest
  id_index <- match(id_name, names(Ttest))
  
  # APPLY FUNCTION
  out<-Ttest[, Tstar.internal(.SD, .BY, conf), by = Ttest[[id_index]]]
  
  return(out)
}

data.TtoStar <- function(Tdata, id_name, conf){
  # ALLOWS USER TO GET data.TstarforId WITHOUT RECALCULATING
  # id_name EQUIVALENT TO THAT OF THE FUNCTION THAT CALC'D TDATA
  Ttest <- Tdata
  Tstar.internal <- function(data, index, conf){
    # UNPACK
    if(is.list(index)) index <- index[[1]]
    
    info <- data[, .SD, .SDcols = c("DateTime", "Tval","Time")]
    dates <- info$DateTime; Ts <- info$Tval; time <- info$Time
    #GET OUT
    Tstar <- as.numeric(max(abs(Ts)))
    t.point <- match(Tstar, abs(Ts))
    # FIND ID
    out <- data.table(DateTime = dates[t.point], Time = time[t.point], Tstar = Tstar)
    out[, paste0(id_name):= index]
    
    if(missing(conf)){
      return(out)
    }
    
    # IF CONF IS ABOVE ONE - THEN WE ASSUME IT IS A THRESHOLD INSTEAD
    mode <- any(conf >= 1)
    if(!mode){
      # FIT RHO
      rho <- est.rho.MLE(Ts)
      rhom <- length(Ts)
      rhorho <- rho
      
      # THE POSSIBILITY OF MULTIPLE CONFIDENCES (can become messy real quick)
      if(length(conf) > 1){
        for(i in 1:length(conf)){
          quantile <- est.z_quantile(rhom, rhorho, conf[i])$qZm
          out[, paste0("q_",conf[i]):= quantile]
          out[, paste0("db_",conf[i]) := Tstar >= quantile]
          
        }
      }
      else{
        quantile <- est.z_quantile(rhom, rhorho, conf)$qZm
        out[, paste0("q_",conf):= quantile]
        out[, paste0("db_",conf) := Tstar >= quantile]
      }
    }
    else{
      # THE POSSIBILITY OF MULTIPLE CONFIDENCES (can become messy real quick)
      if(length(conf) > 1){
        for(i in 1:length(conf)){
          out[, paste0("threshold_",conf[i]):= conf[i]]
          out[, paste0("db_",conf[i]) := Tstar >= conf[i]]
          
        }
      }
      else{
        out[, paste0("threshold__",conf):= conf]
        out[, paste0("db_",conf) := Tstar >= conf]
      }
    }
    return(out)
  }
  # IDENTIFY ID COLUMN IN Ttest
  id_index <- match(id_name, names(Ttest))
  
  # APPLY FUNCTION
  out<-Ttest[, Tstar.internal(.SD, .BY, conf), by = Ttest[[id_index]]]
  
  return(out)
}

# PLOT
data.plot_db<-function(data, burst_time, window = 20*60, hd, hv, lag = 10, blue = F){
  # NEEDS SCALE/OFFSET PARAMETERS FOR AESTETICS
  dt <- data
  # FIND WINDOW
  start<- burst_time-window
  end <- burst_time+window
  start_index <- data.date_To_tindex(start, dt$DateTime)
  end_index <- data.date_To_tindex(end, dt$DateTime)
  
  # EXTRACT DATA
  focus <- dt[start_index:end_index, ]
  
  price <- exp(focus$logPrice)
  Y <- diff(dt$logPrice) # we need all data for estimation
  time <- dt$Time
  DateTime <- focus$DateTime
  
  # ESTIMATE T
  input <- list(time = time, Y = Y)
  t.index <- start_index:end_index       # check this
  
  t.index <- unique(t.index)
  
  # <ESTIMATION>
  mu<-est.mu.next.cpp(input, hd = hd, t.index = t.index)
  sig<-est.sigma.next.cpp(input, hv = hv, t.index = t.index, lag = lag)
  
  # RESCALING
  mu_re <- est.rescale.mu.vec(mu$mu, mu$time, t_beginning = time[1], h_mu = hd)
  sig_re <- est.rescale.sigma.vec(sig$sig, sig$time, t_beginning = time[1], h_sigma = hv)
  
  Ttest <- sqrt(hd)*mu_re/sqrt(sig_re)
  
  # PLOTTERIA
  require(ggplot2)
  plotdata <- data.frame(Date = DateTime, Price = price, T.statistic = Ttest) # mange punkter uden ændring
  
  offset <- mean(price)
  scale <- (max(Ttest)-min(Ttest))/(max(price)-min(price)) # check this later
  
  if(blue){
    p <- ggplot(data = plotdata, aes(x = Date)) +
      geom_area(aes(y = T.statistic)) + #rescale this
      geom_line(aes(y = (Price-offset)*scale, colour = "Price"), color = "dodgerblue3", size = 1) + # size might change
      scale_y_continuous(sec.axis = sec_axis(trans = ~./scale+offset, name = "Price"))
  }
  else{
    p <- ggplot(data = plotdata, aes(x = Date)) +
      geom_area(aes(y = T.statistic)) + #rescale this
      geom_line(aes(y = (Price-offset)*scale, colour = "Price"), size = 1) + # size might change
      scale_y_continuous(sec.axis = sec_axis(trans = ~./scale+offset, name = "Price"))
  }
  
  p
  return(p)
}

data.floor_date <- function(time, unit = "sec"){
  # floors the date time
  # unit can be one of the time units
  decimals <- time - round(time, paste0(unit))
  return(time-decimals)
}

# Knaldhamrende fast - 4200x faster than naive implementation
data.date_To_tindex<-function(dates, t_dates){
  # ASSUMES SORTED
  date <- t_dates
  sorter <- data.table(date, val = date)
  #setattr(sorter, "sorted", "date")
  sorter <- setkey(sorter)
  table <- sorter[J(dates), roll = "nearest"]
  out <- match( table$val, t_dates)
  return(out)
}


# TRADES DISTRIBUTION
data.trades_dist <- function(data, bucketlength, bucketstart, bucketend){
  # RETURNS GGPLOTTABLE DATA
  require(ggplot2)
  data <- data.dayID(data)
  ndays <- length(unique(data$day))
  # Remove date
  t <- strftime(data$DateTime, format = "%H:%M:%S", tz = "UTC")
  
  # reformat as posixct - this forces all to same day (today)
  ntimes <- as.POSIXct(t, format = "%H:%M:%S", tz = "UTC")
  # SORT
  ntimes <- sort(ntimes)
  
  if(missing(bucketstart) || missing(bucketend)){
    # Setup buckets
    buckets <- seq(data.floor_date(ntimes[1], unit = "day"), data.floor_date(ntimes[1], unit = "day")+3600*24, by = bucketlength)[-1]
  }
  else{
    # MARK TERRITORY
    dates<-as.Date(ntimes)[1]
    
    start <- as.POSIXct(paste0(dates, " ", bucketstart), tz = "UTC")
    end <- as.POSIXct(paste0(dates, " ", bucketend), tz = "UTC")
    
    buckets <- seq(start, end, by = bucketlength)[-1]
  }
  
  bucket.id <- .bincode(ntimes, breaks = c(0,buckets))
  
  dt <- data.table(Time = ntimes, bucket = bucket.id)
  
  N <- dt[, .N, by = bucket]$N/ndays
  
  # PLOT DATA
  return(plot.data <- data.table(Time = buckets, N = N))
}

# EXAMPLES
# data.trades_dist(sp, 300, "09:30:00", "16:00:00")
# data.trades_dist(bit, 300)








# REVERSE (only tested daily)
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
      mu <- est.mu.next.cpp(estdata, prevmu, hd, estdata_index)
      sig <- est.sigma.next.cpp(estdata, prevsig, hv, estdata_index, lag = lag) # DEN SKAL IKKE RETURNE PREVSIG
      testT <- sqrt(hd)*mu$mu/sig$sig
      
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

#test <- reverse_ID(dt, grid_T, grid, 7, X = c(300,600), Y = c(300), id_name = id_name)

# REVERSE MATRIX for easier application
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

data.cont <- function(DT, maxTol){
  #extract longest continous period
  t <- DT$Time
  maxTol <- 120
  dt <- diff(t)
  
  ddt <- diff(c(0,which(dt>=maxTol), length(t)))
  
  # crit <- ddt==max(ddt[ddt!=max(ddt)]) # second longest interval
  crit <- ddt==max(ddt) # longest interval. 
  maxPeriod <- which(crit)
  
  StartEnd <-  c(c(0,which(dt>=maxTol), length(t))[maxPeriod] + 1 , c(0,which(dt>=maxTol), length(t))[maxPeriod+1])
  DT2 <- DT[StartEnd[1]:StartEnd[2],]
  
  return(DT2)
}

data.Changed <- function(Vec){
  # finds Changing values
  dn <- diff(Vec)
  return( c(TRUE, dn != 0))
}
