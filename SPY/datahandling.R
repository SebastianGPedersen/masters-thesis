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

data.xsecID <- function(datatable, bucketLengthInSeconds, id = "id"){
  # datatable should be a data.table
  # id denotes the name of the id column
  dt<-daysInRange
  time <- dt[, "DateTime"]
  #Intraday
  buckets <- seq(time$DateTime[1], time$DateTime[length(time$DateTime)], by = bucketLengthInSeconds)
  dtB2<- .bincode(dt$DateTime, breaks = c(0, buckets, time$DateTime[length(time$DateTime)]))
  dt[, paste0(id) := dtB2]
  gc()
  return(dt)
}

data.test.db<-function(time, dy, hd, hv, kern = kern.leftexp, t.index){
  # data should include time | Y (log returns) | raw (before preaverage)
  # Calculates db test stat by computing mu/sig
  
  data <- list(time = time, Y = dy)
  mu<-est.mu.next(data = data, hd = hd, kern = kern, t.index = t.index)
  sig <- est.sigma.next(data = data, hv=hv, t.index = t.index, kern = kern)
  
  # Calculate T
  Tstat<-teststat(mu, sig, hd, hv, kern)
  
  return(Tstat)
}