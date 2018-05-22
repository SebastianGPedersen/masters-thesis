options(digits = 14, digits.secs = 5)

#### source and get data - slow
cd <- getwd()
setwd(paste0(Sys.getenv("masters-thesis"),"/SPY"))
source("dataFunctions.R")
source("datahandling.R")

setwd(paste0(Sys.getenv("masters-thesis-data"),"/SPY"))
fullData<-readRDS("2014_SPY_Vol_Avg.rds")
setwd(cd)
rm(cd)

#### various data extractions
#Date format: yyyy-mm-dd

# daysInRange<-selectDays(fullData, as.Date("2014-01-01"), endDate = as.Date("2014-01-03"))
daysInRange<-fullData
rm(fullData)


###################### CREATE DAY ID #############################

data.dayID <- function(datatable, id = "id"){
  # datatable should be a data.table.
  # id denotes the name of the id column
  dt <- datatable
  # sort by dateTime
  setkey(dt, DateTime)
  # Identify unique days
  unDates <- unique(dt$Date)
  # Match with ID number
  dtB<-match(dt$Date, unDates)
  # Add ID column
  dt[, paste0(id) := dtB]
  gc()
  return(dt)
}

if(F){
  dt<-daysInRange
  # sort by dateTime
  setkey(dt, DateTime)
  
  gg <- "uld"
  
  # Identify IDs
  unDates <- unique(dt$Date)
  dtB<-match(dt$Date, unDates)
  dt[, paste0(gg) := dtB]
  data <- dt[, -c("DateTime", "Date")]
  rm(dtB)
  gc()
}

#################### CREATE X-SECOND BIN ##########################

data.xsecID <- function(datatable, bucketLengthInSeconds, id = "id"){
  dt<-daysInRange
  time <- dt[, "DateTime"]
  #Intraday
  buckets <- seq(time$DateTime[1], time$DateTime[length(time$DateTime)], by = bucketLengthInSeconds)
  dtB2<- .bincode(dt$DateTime, breaks = c(0, buckets, time$DateTime[length(time$DateTime)]))
  dt[, paste0(id) := dtB2]
  gc()
  return(dt)
}

if(F){
  dt<-daysInRange
  time <- dt[, "DateTime"]
  #Intraday
  buckets <- seq(floor_date(time$DateTime[1]), time$DateTime[length(time$DateTime)], by = "hour")
  dtB2<- .bincode(dt$DateTime, breaks = c(buckets, time$DateTime[length(time$DateTime)]))
  dt[, "hour" := dtB2]
  gc()
  # food4thought
  #N <- dt[, .N, by = id]
  #plot(N$N)
}

#################### APPLY FCT PER ID ##########################

dt <- daysInRange
dt <- data.dayID(dt, "day")

dy <- dt[, lapply(.SD, diff), by = day, .SDcols = c("logPrice")]$logPrice
time <- dt[, "Time"]$Time
max(diff(time))


#################### THOUGHTS ####################

# REAL DATA (non-bitcoin) are broken by nights - therefore any estimation should only be done on observations of same day

# BITCOIN DATA there are no nights - estimation should be done on all previous observations - (costly)

################### /THOUGHTS ####################


floor_date <- function(time, unit = "sec"){
  # floors the date time
  # unit can be one of the time units
  decimals <- time - round(time, paste0(unit))
  return(time-decimals)
}

# START # START # # START # START # # START # START # # START # START # # START # START # 
# START # START # # START # START # # START # START # # START # START # # START # START # 
# START # START # # START # START # # START # START # # START # START # # START # START # 

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

# APPLY ON DAY #

data<-data.getFull()
data<-data.dayID(data, id = "day")

# DAY FUNCTION START
data.TforId <- function(data, id_name, id_number, hd, t.freq, lag, offset){
  # Estimates T test for every t.freq (unit: seconds) within a certain id
  # hd and lag are parameters for estimation
  # offset is how many observations are skipped

  # IDENTIFY COLUMN INDEX
  if(id_name %in% names(data)){
    id_index <- match(id_name, names(data))
  }
  else{
    stop("id_name not found in data")
  }
  
  # EXTRACT RELEVANT ID
  .id <- id_number
  dt<-data[data[[id_index]] == .id, ]
  
  #print(id_number)
  
  # unpack data
  dy <- dt[, lapply(.SD, diff), .SDcols = "logPrice"]$logPrice
  times <- dt[, .SD, .SDcols = c("DateTime","Time")]
  time <- times$Time; Date <- times$DateTime
  
  .ids <- dt[, .SD, .SDcols = c(id_name)][[1]]
  
  # DETERMINE INPUTS
  input <- list(time = time, Y = dy)
  
  # SETUP T-INDEX #
  dates<-seq(data.floor_date(Date[1]), Date[length(Date)], by = t.freq)
  t.index<-data.date_To_tindex(dates, t_dates = Date)
  
  # OFFSET
  t.index <- t.index[(1+offset):length(t.index)] # BUT THIS ONLY SKIPS ON DAY 1
  
  # <ESTIMATION>
  #mu<-est.mu.next(input, hd = hd, t.index = t.index)$mu
  #sig<-est.sigma.next(input, hv = hd, t.index = t.index, lag = lag)$sig # t.index does not work because last is included - THINK ABOUT THIS
  #Ttest<-sqrt(hd/kern.leftexp$ksq)*mu/sqrt(sig)
  # </ESTIMATION>
  
  mu <- 1:length(t.index)
  sig <- 1:length(t.index)*2        # testeria
  Ttest <- 1:length(t.index)*10
  
  # FORMULATE OUTPUT
  Nout<-length(t.index)
  out <- data.table(DateTime = Date[t.index], Time = time[t.index], Mu = mu, Sigma = sig, Tval = Ttest)
  out[, paste0(id_name) := .ids[t.index]]
}

res <- data.TforId(data = data, id_name = "day", id_number = 2, hd = 1, t.freq = 60, lag = 10, offset = 10)

#require(microbenchmark)
#microbenchmark(data.TforId(data = data, id_name = "day", id = 1, hd = 1, t.freq = 60, lag = 10, offset = 10), times = 30)<- 1.3s pr day

if(F){
  # get diff
  .id <- id
  dt<-data[day == .id, ]
  #dt<-data[day <= .id, ]
  
  dy <- dt[, lapply(.SD, diff), .SDcols = "logPrice"]$logPrice
  times <- dt[, .SD, .SDcols = c("DateTime","Time")]
  time <- times$Time; Date <- times$DateTime
  
  .ids <- dt[, .SD, .SDcols = c("day")]$day
  
  # DETERMINE INPUTS
  input <- list(time = time, Y = dy)
  
  # SETUP T-INDEX #
  dates<-seq(data.floor_date(Date[1]), Date[length(Date)], by = t.freq)
  t.index<-data.date_To_tindex(dates, t_dates = Date)
  
  # SKIP START
  t.index <- t.index[(1+offset):length(t.index)] # BUT THIS ONLY SKIPS ON DAY 1
  
  # <ESTIMATION>
  #mu<-est.mu.next(input, hd = hd, t.index = t.index)$mu
  #sig<-est.sigma.next(input, hv = hd, t.index = t.index, lag = lag)$sig # t.index does not work because last is included - THINK ABOUT THIS
  #Ttest<-sqrt(hd/kern.leftexp$ksq)*mu/sqrt(sig)
  # </ESTIMATION>
  
  mu <- 1:length(t.index)
  sig <- 1:length(t.index)*2        # testeria
  Ttest <- 1:length(t.index)*10
  
  # FORMULATE OUTPUT
  Nout<-length(t.index)
  out <- data.table(DateTime = Date[t.index], Time = time[t.index], Mu = mu, Sigma = sig, Ttest = Ttest, day = .ids[t.index])
}

# MAKE ABOVE A FUNCTION SUCH THAT IT CAN BE APPLIED OVER EVERY ' by = id '

subdata <- data[day <= 5 ,]

res <- data.TforId(data = subdata, id_name = "day", id_number = 1, hd = 1, t.freq = 60, lag = 10, offset = 10)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# CALCULATE TSTAR ON OUTPUT # (should be new function, mmmhyes?)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data.TstarforId <- function(data, id_name, id_number, hd, t.freq, lag, offset, conf){
  # WRAPPER AROUND data.TforId THAT RETURNS THE TSTAR VALUE INSTEAD
  
  Ttest <- data.TforId(data = data, id_name = id_name, id_number = id_number, hd = hd, t.freq = t.freq, lag = lag, offset = offset)
  # UNPACK
  info <- Ttest[, .SD, .SDcols = c("DateTime", "Tval","Time")]
  dates <- info$DateTime; Ts <- info$Tval; time <- info$Time
  #GET OUT
  Tstar <- as.numeric(max(abs(Ts)))
  index <- match(Tstar, Ts)
  # FIND ID
  id_index <- match(id_name, names(data))
  ids <- Ttest[, .SD, .SDcols = c(id_name)][[1]]
  out <- data.table(DateTime = dates[index], Time = time[index], Tstar = Tstar)
  out[, paste0(id_name):= ids[index]]
  
  if(missing(conf)){
    return(out)
  }
  
  # FIT RHO
  rho <- est.rho(Ts)
  rhom <- rho$m
  rhorho <- rho$rho
  
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
  
  return(out)
}

data.TtoStar <- function(Tdata, id_name, conf){
  # ALLOWS USER TO GET data.TstarforId WITHOUT RECALCULATING
  # id_name EQUIVALENT TO THAT OF THE FUNCTION THAT CALC'D TDATA
  Ttest <- Tdata
  # UNPACK
  info <- Ttest[, .SD, .SDcols = c("DateTime", "Tval","Time")]
  dates <- info$DateTime; Ts <- info$Tval; time <- info$Time
  #GET OUT
  Tstar <- as.numeric(max(abs(Ts)))
  index <- match(Tstar, Ts)
  out <- data.table(DateTime = dates[index], Time = time[index], Tstar = Tstar)
  # FIND ID
  id_index <- match(id_name, names(data))
  ids <- Ttest[, .SD, .SDcols = c(id_name)][[1]]
  out[, paste0(id_name):= ids[index]]
  
  if(missing(conf)){
    return(out)
  }
  
  # FIT RHO
  rho <- est.rho(Ts)
  rhom <- rho$m
  rhorho <- rho$rho
  
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
  return(out)
}

res

dild<-data.TstarforId(subdata, id_name = "day", id_number = 1, hd = 1, t.freq = 60, lag = 10, offset = 10, conf = 0.95)
dild

tst <- data.TtoStar(res, "day", conf = c(0.95, 0.995))
dild
tst

#dt<-data.xsecID(dt, 3600, id = "egg") # We have to name it something that no one else will
#dt<-dt[egg == id, ]

#subdata[, data.TstarforId(subdata, id_name = "day", id_number = day, hd = 1, t.freq = 60, lag = 10, offset = 10), by =  .I] NO :(

#######################################################################################################

# ALLOW ID INTERVAL FROM A TO B
# option to offset for every new id? (this is basically a forloop of the old with table merge at end) #

#######################################################################################################


#id_name <- "day"; hd <- 1; t.freq <- 60; lag <- 10; offset <- 100; offset_perId <- T
#id_number <- 2:3

# THIS ONE IS BAD
data.TforId <- function(data, id_name, id_number, hd, t.freq, lag, offset){
  # Estimates T test for every t.freq (unit: seconds) within a certain id
  # hd and lag are parameters for estimation
  # offset as vector, where each element is how many observations are skipped for that id
  # If offset_perId is false, then offset is only applied on first
  # If offset_perId is true, but offset is only a single value, then this value will be used for every id
  # Bandwidth, freq, lag and offset can all be vectors
  
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
  
  # IDENTIFY COLUMN INDEX
  if(id_name %in% names(data)){
    id_index <- match(id_name, names(data))
  }
  else{
    stop("id_name not found in data")
  }
  
  # DEFINE INTERNAL FUNCTION FOR EASY LOOP
  TforId.internal<-function(data, id_name, id_number, hd, t.freq, lag, offset){
    # EXTRACT RELEVANT ID
    dt<-data[data[[id_index]] == id_number, ]
    
    
    # unpack data
    dy <- dt[, lapply(.SD, diff), .SDcols = "logPrice"]$logPrice
    times <- dt[, .SD, .SDcols = c("DateTime","Time")]
    time <- times$Time; Date <- times$DateTime
    
    #ids <- dt[, .SD, .SDcols = c(id_name)][[1]]
    ids <- dt[[id_name]]
    
    # DETERMINE INPUTS
    input <- list(time = time, Y = dy)
    
    # SETUP T-INDEX #
    dates<-seq(data.floor_date(Date[1]), Date[length(Date)], by = t.freq)
    t.index<-data.date_To_tindex(dates, t_dates = Date)
    
    # OFFSET
    t.index <- t.index[(1+offset):length(t.index)] # BUT THIS ONLY SKIPS ON DAY 1
    
    # <ESTIMATION>
    #mu<-est.mu.next(input, hd = hd, t.index = t.index)$mu
    #sig<-est.sigma.next(input, hv = hd, t.index = t.index, lag = lag)$sig # t.index does not work because last is included - THINK ABOUT THIS
    #Ttest<-sqrt(hd/kern.leftexp$ksq)*mu/sqrt(sig)
    # </ESTIMATION>
    
    mu <- 1:length(t.index)
    sig <- 1:length(t.index)*2        # testeria
    Ttest <- 1:length(t.index)*10
    
    # FORMULATE OUTPUT
    Nout<-length(t.index)
    out <- data.table(DateTime = Date[t.index], Time = time[t.index], Mu = mu, Sigma = sig, Tval = Ttest)
    out[, paste0(id_name) := ids[t.index]]
  }
  
  out <- TforId.internal(data, id_name = id_name, id_number = id_number[1],
                         hd = hd[1], t.freq = t.freq[1], lag = lag[1], offset = offset[1])
  # HERE WE LOOP OVER THE IDs
  if(length(id_number > 1)){
    for(i in 2:length(id_number)){
      run <- TforId.internal(data, id_name = id_name, id_number = id_number[i],
                             hd = hd[i], t.freq = t.freq[i], lag = lag[i], offset = offset[i])
      # BIND TOGETHER
      out <- rbind(out, run)
    }
  }
  return(out)
}

# ATTEMPT AT BETTER DATATABLE-ID HANDLING

# OLD - SLIGHTLY IMPROVED IN THE DATAHANDLING.R
data.TforId <- function(data, id_name, id_number, hd, t.freq, lag, offset, offset_perId = T){
  # Estimates T test for every t.freq (unit: seconds) within a certain id
  # hd and lag are parameters for estimation
  # offset as vector, where each element is how many observations are skipped for that id
  # If offset_perId is false, then offset is only applied on first
  # If offset_perId is true, but offset is only a single value, then this value will be used for every id
  # Bandwidth, freq, lag and offset can all be vectors
  
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
  
  # IDENTIFY COLUMN INDEX
  if(id_name %in% names(data)){
    id_index <- match(id_name, names(data))
  }
  else{
    stop("id_name not found in data")
  }
  
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
    
    # OFFSET
    t.index <- t.index[(1+offset):length(t.index)] # BUT THIS ONLY SKIPS ON DAY 1
    
    # <ESTIMATION>
    #mu<-est.mu.next(input, hd = hd, t.index = t.index)$mu
    #sig<-est.sigma.next(input, hv = hd, t.index = t.index, lag = lag)$sig # t.index does not work because last is included - THINK ABOUT THIS
    #Ttest<-sqrt(hd/kern.leftexp$ksq)*mu/sqrt(sig)
    # </ESTIMATION>
    
    mu <- 1:length(t.index)
    sig <- 1:length(t.index)*2        # testeria
    Ttest <- 1:length(t.index)*10
    
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

res <- data.TforId(data = data, id_name = "day", id_number = 2:3, hd = 1, t.freq = c(30,60), lag = 10, offset = 10, offset_perId = T)

data.TstarforId <- function(data, id_name, id_number, hd, t.freq, lag, offset, conf){
  # WRAPPER AROUND data.TforId THAT RETURNS THE TSTAR VALUE INSTEAD
  
  
  Ttest <- data.TforId(data = data, id_name = id_name, id_number = id_number, hd = hd, t.freq = t.freq, lag = lag, offset = offset)
  
  Tstar.internal <- function(data, index, conf){
    # UNPACK
    if(is.list(index)) index <- index[[1]]
    
    info <- Ttest[, .SD, .SDcols = c("DateTime", "Tval","Time")]
    dates <- info$DateTime; Ts <- info$Tval; time <- info$Time
    #GET OUT
    Tstar <- as.numeric(max(abs(Ts)))
    t.point <- match(Tstar, Ts)
    # FIND ID
    out <- data.table(DateTime = dates[t.point], Time = time[t.point], Tstar = Tstar)
    out[, paste0(id_name):= index]
    
    if(missing(conf)){
      return(out)
    }
    
    # FIT RHO
    rho <- est.rho(Ts)
    rhom <- rho$m
    rhorho <- rho$rho
    
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
    return(out)
  }
  # IDENTIFY ID COLUMN IN Ttest
  id_index <- match(id_name, names(Ttest))
  
  # APPLY FUNCTION
  out<-Ttest[, Tstar.internal(.SD, .BY, conf), by = Ttest[[id_index]]]
  
  return(out)
}

dild<-data.TstarforId(data, id_name = "day", id_number = 1:3, hd = 1, t.freq = 60, lag = 10, offset = 10, conf = 0.95)
dild

data.TtoStar <- function(Tdata, id_name, conf){
  # ALLOWS USER TO GET data.TstarforId WITHOUT RECALCULATING
  # id_name EQUIVALENT TO THAT OF THE FUNCTION THAT CALC'D TDATA
  Ttest <- Tdata
  Tstar.internal <- function(data, index, conf){
    # UNPACK
    if(is.list(index)) index <- index[[1]]
    
    info <- Ttest[, .SD, .SDcols = c("DateTime", "Tval","Time")]
    dates <- info$DateTime; Ts <- info$Tval; time <- info$Time
    #GET OUT
    Tstar <- as.numeric(max(abs(Ts)))
    t.point <- match(Tstar, Ts)
    # FIND ID
    out <- data.table(DateTime = dates[t.point], Time = time[t.point], Tstar = Tstar)
    out[, paste0(id_name):= index]
    
    if(missing(conf)){
      return(out)
    }
    
    # FIT RHO
    rho <- est.rho(Ts)
    rhom <- rho$m
    rhorho <- rho$rho
    
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
    return(out)
  }
  # IDENTIFY ID COLUMN IN Ttest
  id_index <- match(id_name, names(Ttest))
  
  # APPLY FUNCTION
  out<-Ttest[, Tstar.internal(.SD, .BY, conf), by = Ttest[[id_index]]]
  
  return(out)
}

tst <- data.TtoStar(res, "day", conf = c(0.95, 0.995))
dild
tst
# TEST THE SEVERAL IDs FUNCTION BY MANUALLY DOING SHIT
# MAKE TSTAR WORK WITH SEVERAL ID NUMBERS! (SHOULDNT TAKE LONG)

# PRO-LEAGUE: MAKE FUNCTION SUCH THAT IT CAN TaKE THE NEEDED FROM PREVIOUS LAG (depends on hd)
#             SUCH THAT OFFSET IS NOT NEEDED IN EVERY ID (THINK HOURLY OR BITCOIN)

# ALLOW FOR HD/HV
