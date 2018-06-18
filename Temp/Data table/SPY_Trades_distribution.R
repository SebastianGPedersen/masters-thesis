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
# rm(fullData)

###################### New stuff #############################
dt<-daysInRange[,-c("Time","logPrice")]
rm(daysInRange)
rm(fullData)
setkey(dt, DateTime)

#create 5min buckets
allDays <- as.POSIXct.Date(as.Date(unique(dt$Date))) 

# allDays <- as.POSIXct.Date(seq(from = as.Date(min(dt$Date)), to =as.Date(max(dt$Date)), by = 1)) #Could be replaced with something like as.POSIXct.Date(unique(dt$Date)) 
buckets2 <- seq(from = 34200+300, to = 57600, by = 300)
buckets <- .POSIXct(character(length(allDays)*length(buckets2)), tz = "UTC") #rep(NA, length(allDays)*length(buckets2))

for(i in seq_along(allDays)){
  buckets[(((i-1)*length(buckets2)):(i*length(buckets2)))[-1]] <- allDays[i]+buckets2+0.0001
}

#Slowly extract share of oberservations in each bucket
setkey(dt, DateTime)

X <- dt[ , list(DailyCount = .N) , by = "Date" ]
# dt2 <- dt[ , list(Date, Bucket = which.max(DateTime <= buckets) ) , by = "DateTime" ]
dt2 <- dt

dt2B<-.bincode(dt$DateTime, breaks = c(0, buckets))

# all.equal(dt2B, dt2B_old)

dt2[, "Bucket" := dt2B]
setkey(dt2, Date)

#Better
dttemp <- dt2[ , list(Count = .N, Date) , by = list(Bucket) ]
setkey(dttemp, Date)
dttemp <- unique(dttemp, by = "Bucket")
dttemp <- dttemp[X, list(Date, Bucket, Count, DailyCount, Share = Count/DailyCount, IntraDayBucket = Bucket %% length(buckets2))]
dttemp$IntraDayBucket[dttemp$IntraDayBucket==0] <- length(buckets2)
dttemp

avgShareIntradayBucket <- dttemp[, list(Share = mean(Share)), by = IntraDayBucket]
bpddt <- dttemp[, list(bpd = .N), by = Date]
plot(bpddt$bpd) #number of buckets per day

plot(avgShareIntradayBucket)  # distribution of trades: Average share of trades within each 5 min interval
plot(dttemp[,mean(Share*DailyCount), by = IntraDayBucket]$V1) # Same as above, but the values rather than the ratio

trades_dist <- dttemp[,mean(Share*DailyCount), by = IntraDayBucket]$V1

plot(trades_dist)


# SAVE
setwd(Sys.getenv("masters-thesis"))
# PREP FOR EXPORT
time <- buckets[1:length(trades_dist)]

trades <- data.frame(time = time ,trades = trades_dist)
save(trades, file = "Simulation/trades.RDa")



#
#
#
# ## Invesigate stuff
# require(lubridate)
# 
# 
# bpddt$Date[which(bpddt$bpd<=70)]
# dt[dt$Date=="2014-03-31"]
# plot(avgShareIntradayBucket)
# dttemp[dttemp$Date==unique(dttemp$Date)[4]]
# buckets[6930]
# dt[dt$Date=="2014-03-31"]$DateTime[length(dt[dt$Date=="2014-03-31"]$DateTime)]<buckets[6930]
