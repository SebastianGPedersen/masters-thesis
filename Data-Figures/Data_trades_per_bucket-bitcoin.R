# BITCOIN DATA TESTING
setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/rescaling.R")
source("kernels/kernels.R")
source("spy/dataFunctions.R")
source("spy/datahandling.R")

require(data.table)

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

bitfinex <- data.getbitdata("bitfinex")
bitfinex <- data.trades_dist(bitfinex, 300)

bitmex <- data.getbitdata("bitmex")
bitmex <- data.trades_dist(bitmex, 300)

kraken <- data.getbitdata("kraken")
kraken <- data.trades_dist(kraken, 300)

plot.data <- data.table(Time = bitfinex$Time, N_bit = bitfinex$N, N_mex <- bitmex$N, N_kraken <- kraken$N)

buckets <- plot.data$Time
timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])

labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")
# PLOT
g1 <- ggplot(data = plot.data) +
  geom_point(aes(x = Time, y = N_bit, group = 1, col = "BitFinex")) +
  geom_point(aes(x = Time, y = N_mex, group = 1, col = "BitMex")) +
  geom_point(aes(x = Time, y = N_kraken, group = 1, col = "Kraken")) +
  ylab("Trades per 5min bucket") +
  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)

g1
