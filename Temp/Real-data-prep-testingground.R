options(digits = 14, digits.secs = 5)

#### source and get data - slow
cd <- getwd()
setwd(paste0(Sys.getenv("masters-thesis"),"/SPY"))
source("dataFunctions.R")

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
  buckets <- seq(time$DateTime[1], time$DateTime[length(time$DateTime)], by = "hour")
  dtB2<- .bincode(dt$DateTime, breaks = c(0, buckets, time$DateTime[length(time$DateTime)]))
  dt[, "id" := dtB2]
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


