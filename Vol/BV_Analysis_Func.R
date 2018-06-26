cd_UNIQUE_NAME_BV_ANALYSIS__ <- getwd()

setwd(paste0(Sys.getenv("masters-thesis"), "\\Vol"))
source("vol_estimators.R")
source("FlexibleFourierForm_Func.R")

setwd(cd_UNIQUE_NAME_BV_ANALYSIS__)
rm(cd_UNIQUE_NAME_BV_ANALYSIS__)


BV.get_SPY_data <- function(maxDate = "2014-12-31"){
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


#########  INIT ##############
dt<-daysInRange[,-c("Time")]
rm(daysInRange)


setkey(dt, DateTime)
dt<-dt[dt$Date != "2014-10-30"]

dt<-dt[dt$Date<=maxDate]

return(dt)
}

BV.data_deseason_BV_Func <- function(dt, bucketLengthInMinutes){
#########################   Bucket stuff ####################
buckets <- vol.est.DataIntradayBucket(DT = dt, m = bucketLengthInMinutes)
dtB2<- .bincode(dt$DateTime, breaks = c(0, buckets))
dt[, "id" := dtB2]
rm(dtB2)
# rm(dt)
gc()
print("IV estimation reached (SLOW)")
###################### IV Intraday stuff #############################
resListIntraday <- vol.est.IVestPathwise(dt)
print("IV estimation completed")
###################### Vol (from IV) stuff ############################
########## BVS
bvSDT <- resListIntraday$IVest$bvSDT
bvSDT
# YEARLY parametrization
bvSDT[, Vol:= sqrt(bvS1 * bucketLengthInMinutes/(60*24*252)),]

nIntradayBuckets <- 390/bucketLengthInMinutes
bvSDT[, IntraDayBucket := id %% nIntradayBuckets]
bvSDT$IntraDayBucket[bvSDT$IntraDayBucket==0] <- nIntradayBuckets

### FFF way
bvSDTfff <- copy(bvSDT)
bvSDTfff[, c("DayBucket", "LogVol") := list(as.factor(as.numeric(as.factor(as.Date(as.character(buckets[bvSDTfff$id]))))), log(Vol))]

res <- FFF.estOptimalP(DT = bvSDTfff, maxP =  nIntradayBuckets/2 - 1, dataCol = "LogVol", IntradayCol = "IntraDayBucket", dayCol = "DayBucket", dailyInteraction = F, dailyParameter = F)

##### Alternate simple way - better for Jensen's inequality 
### Taking log first to bring in line with fff
fit <- lm(bvSDTfff$LogVol ~ 1 + as.factor(bvSDTfff$IntraDayBucket))
simpleCorrection <- FFF.singlePeriodFit(fit = fit,DT = bvSDTfff, IntradayCol = "IntraDayBucket")

correctionNice <- FFF.singlePeriodFit(res$fit, res$DT, IntradayCol = "IntraDayBucket")

########################### ACF  #########################
bvSDTfff[, sCorrect := simpleCorrection[IntraDayBucket],]
bvSDTfff[, LogVolsCorrect := LogVol - sCorrect,]
bvSDTfff[, nCorrect := correctionNice[IntraDayBucket],]
bvSDTfff[, LogVolnCorrect := LogVol - nCorrect,]
IntradayPerDay <- (60/bucketLengthInMinutes)*24
bvSDTfff[, lagInd:= (as.numeric(DayBucket)-1)*IntradayPerDay + (IntraDayBucket)]

return(list(bvSDTfff=bvSDTfff, optimalP=res$optimalP, buckets=buckets, countDT = resListIntraday$countDT))

}