# BITCOIN DATA PREPARATION #

setwd(Sys.getenv("masters-thesis"))
source("spy/dataFunctions.R")
source("spy/datahandling.R")

require(data.table)

options(digits = 16, digits.secs = 6)

# REORGANIZE DATA AS RDS
setwd(paste0(Sys.getenv("masters-thesis-data"),"/Bitcoin"))

data.prep<-function(data, name){

  
  # NAME FILE
  drop<-paste0(name,".rds")
  # LOGPRICE
  logprice <-log(rowMeans(data.table(ask = data$ask_price, bid = data$bid_price)))
  
  # DATE-TIME
  datetime <- as.POSIXct(data$time_stamp, tz = "UTC")
  
  # SORT DATA
  sortable <- data.table(Datetime = datetime, logPrice = logprice)
  setkey(sortable, sortby = Datetime)
  datetime <- sortable$Datetime
  
  # NON UNIQUE HANDLING
  data.equidistant<-function(datetime){
    # Distributes observations across the previous second
    # such that (2,2) -> (1.5, 2)
    breaks <- c(0,unique(datetime))
    timepoint <- .bincode(datetime, breaks = breaks)
    
    # force double by using miliseconds                                 # holy shit this is mindblowing
    dt <- data.table(Datetime = datetime, Timepoint = timepoint)
    
    equidistant <- function(data){
      mid <- data
      if(length(mid$Datetime)==1){
        return() 
      }
      new<-seq(from = mid$Datetime[1]-1, to = mid$Datetime[1], length.out = length(mid$Datetime)+1)[-1]
      return(as.POSIXct(new))
    }
    
    plurals<-dt[, equidistant(.SD), by = Timepoint]
    
    plurals <- data.table(Timepoint = plurals$Timepoint, Datetime = plurals$V1)
    
    all<-unique(rbind(plurals, unique(dt)))
    setkey(all, sortby = Datetime)
    
    return(all$Datetime)
  }
  datetime <- data.equidistant(datetime = datetime) # takes ~200 seconds - ggwp
  
  # DATE
  date <- as.Date(datetime)
  
  # TIME
  time <- as.numeric(datetime)
  
  # COMBINE
  dt <- data.table(Date = date, DateTime = datetime, Time = time, logPrice = logprice)
  
  # checks
  if(length(dt$DateTime)-length(unique(dt$DateTime)) != 0) stop("Unique handling went wrong - Datetime")
  if( length(dt$Time)-length(unique(dt$Time)) != 0) stop("Unique handling went wrong - Time")
  
  # DELIVERY
  saveRDS(dt, file = drop)
}

# --- BITFINEX ---
load("bitfinex.rda")

data.prep(bitfinex, "bitfinex")

# --- BITMEX ---
load("bitMEX.rda")

data.prep(bitMEX, "bitmex")

# --- KRAKEN ---
load("Kraken.rda")

data.prep(Kraken, "kraken")

# CHECKS
data.getdata <- function(name){
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

bitfigner <- data.getdata("bitfinex")

bitmexico <- data.getdata("bitmex")

kraken <- data.getdata("kraken")

