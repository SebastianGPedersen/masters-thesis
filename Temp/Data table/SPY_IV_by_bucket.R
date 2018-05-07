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


###################### Buckets Day stuff #############################
dt<-daysInRange[,-c("Time")]
rm(daysInRange)


setkey(dt, DateTime)
# dt<-dt[dt$Date != "2014-10-30"]

#Full day
unDates <- unique(dt$Date)
dtB<-match(dt$Date, unDates)
dt[, "id" := dtB]
rm(dtB)
gc()
###################### IV Day stuff #############################

# resList <- vol.est.IVestPathwise(dt)
# resList <- readRDS("FullDayresList.rds")
# saveRDS(resList, "FullDayresList.rds")

resList$IVest

mIvest<-c(mean(resList$IVest$rvDT$rv),
          mean(resList$IVest$rvSDT$rvS),
          mean(resList$IVest$bvDT$bv),
          mean(resList$IVest$bvSDT$bvS))
names(mIvest) <- c("rv", "bv", "rvS", "bvS")
mIvest

plot(as.Date(unDates),resList$IVest$rvDT$rv, xlab = "Date", ylab = "", main = "Daily IV estimates",  type = "l", col = "black")
lines(as.Date(unDates),resList$IVest$rvSDT$rvS, col = "red")
lines(as.Date(unDates),resList$IVest$bvDT$bv, col = "blue")
lines(as.Date(unDates),resList$IVest$bvSDT$bvS, col = "green")
###################### Omega2 #############################
mean(unlist(resList$omega2est))
mean(sqrt(unlist(resList$omega2est)))
###################### Noise ratio ############################

BucketCount <- dt[,.N,by="id"]$N
noiseRatiorvS <- sqrt(BucketCount*resList$omega2est/resList$IVest$rvSDT$rvS)
noiseRatiobvS <- sqrt(BucketCount*resList$omega2est/resList$IVest$bvSDT$bvS)
quantile(noiseRatiobvS$r1)

###################### Vol (from IV) stuff ############################
sqrt(mIvest*(6.5/(252*24))^(-1))
sqrt(0.15*252)

sqrt(252*resList$IVest$rvDT$rv)
sqrt(252*resList$IVest$bvDT$bv)



###################### Buckets Intraday stuff #############################
#Same dt as above


#Intraday
bucketLengthInMinutes <- 5
buckets <- vol.est.DataIntradayBucket(DT = dt, m = bucketLengthInMinutes)
dtB2<- .bincode(dt$DateTime, breaks = c(0, buckets))
dt[, "id" := dtB2]
rm(dtB2)
# rm(dt)
gc()
###################### IV Intraday stuff #############################
# resListIntraday <- vol.est.IVestPathwise(dt)
resListIntraday <- readRDS("5minresListIntraday.rds")
resListIntraday$IVest

mIvest<-c(mean(resListIntraday$IVest$rvDT$rv),
          mean(resListIntraday$IVest$rvSDT$rvS),
          mean(resListIntraday$IVest$bvDT$bv),
          mean(resListIntraday$IVest$bvSDT$bvS))
plot(resListIntraday$IVest$rvDT$rv, xlab = "Date", ylab = "", main = "5-min IV estimates",  type = "l", col = "black")
lines(resListIntraday$IVest$rvSDT$rvS, col = "red")
lines(resListIntraday$IVest$bvDT$bv, col = "blue")
lines(resListIntraday$IVest$bvSDT$bvS, col = "green")

# saveRDS(resListIntraday, "5minresListIntraday.rds")

# plot(bvSDT[IntraDayBucket==16]$Vol/avg5MinVol)
###################### Vol (from IV) stuff ############################
########## BVS
bvSDT <- resListIntraday$IVest$bvSDT
bvSDT
# YEARLY parametrization
bvSDT[, Vol:= sqrt(bvS1 * bucketLengthInMinutes/(60*24*252)),]
plot(bvSDT$id, bvSDT$Vol)
#78 = 5 min
nIntradayBuckets <- 390/bucketLengthInMinutes
bvSDT[, IntraDayBucket := id %% nIntradayBuckets]
bvSDT$IntraDayBucket[bvSDT$IntraDayBucket==0] <- nIntradayBuckets

### Simple way
avgIntradayVol <- bvSDT[, list(AvgIntradayVol = mean(Vol)), by = IntraDayBucket]

#Slighly inaccurate, as there are NOT the same number of obs per day
avgXMinVol <- mean(avgIntradayVol$AvgIntradayVol)
# simpleCorrection <- avgIntradayVol$AvgIntradayVol/avgXMinVol
# simpleCorrection <- avgIntradayVol$AvgIntradayVol
simpleCorrection <- log(avgIntradayVol$AvgIntradayVol)
plot(avgIntradayVol$IntraDayBucket, simpleCorrection, type = "l")

# plot(bvSDT[DayBucket==220]$Vol/avgXMinVol)

### FFF way
bvSDTfff <- copy(bvSDT)
bvSDTfff[, c("DayBucket", "LogVol") := list(as.factor(as.numeric(as.factor(as.Date(as.character(buckets[bvSDTfff$id]))))), log(Vol))]
# bvSDT[ ,"DayBucket":=bvSDTfff[["DayBucket"]]]

# intraDayI <- bvSDTfff$IntraDayBucket
res <- vol.est.optimalP(DT = bvSDTfff, maxP =  nIntradayBuckets/2 - 1, dataCol = "LogVol", IntradayCol = "IntraDayBucket", dayCol = "DayBucket", dailyInteraction = F, dailyParameter = F)
res$CritVec
res$fit$rank
summary(res$fit)
sum(is.na(res$fit$coefficients))
res$optimalP
plot(res$DT$LogVol[1:nIntradayBuckets])
lines(res$fit$fitted.values[1:nIntradayBuckets])
plot(res$fit$fitted.values[1:nIntradayBuckets]-res$fit$coefficients[1])

plot(exp(vol.est.singlePeriodFit(res$fit, res$DT, IntradayCol = "IntraDayBucket")), type = "l")

# 
# 
# res2 <- vol.est.optimalP(DT = bvSDTfff, maxP = nIntradayBuckets/2 - 1 , dataCol = "LogVol", IntradayCol = "IntraDayBucket", dayCol = "DayBucket", dailyInteraction = F, dailyParameter = T)
# summary(res2$fit)
# sum(is.na(res2$fit$coefficients))
# res2$optimalP   
# plot(exp(vol.est.singlePeriodFit(res2$fit, res2$DT, IntradayCol = "IntraDayBucket"))/avgXMinVol, type = "l")
# lines((res2$fit$fitted.values))
# plot((res2$fit$fitted.values-res2$fit$coefficients[1])[1:78])
# res2$fit$coefficients[1:10]
# plot((res2$fit$fitted.values)[1:78])
# 
# 
# res3 <- vol.est.optimalP(DT = bvSDTfff, maxP = res2$optimalP, dataCol = "LogVol", IntradayCol = "IntraDayBucket", dayCol = "DayBucket", dailyInteraction = F, dailyParameter = F, minP=res2$optimalP   )
# plot(exp(vol.est.singlePeriodFit(res3$fit, res3$DT, IntradayCol = "IntraDayBucket")), type = "l")
# summary(res3$fit)
# plot(exp(res3$fit$fitted.values[1:78]-res3$fit$coefficients[1]))
# 
# plot(res3$fit$coefficients[-1])

##### PLOT
plot(avgIntradayVol$IntraDayBucket, simpleCorrection, type = "l")#, ylim = c(0.6,1.8))
correctionNice <- vol.est.singlePeriodFit(res$fit, res$DT, IntradayCol = "IntraDayBucket")
lines(correctionNice, type = "l", col = "red")
# lines(exp(vol.est.singlePeriodFit(res2$fit, res2$DT, IntradayCol = "IntraDayBucket"))/avgXMinVol, type = "l", col = "blue")
# correctionOverfit <- exp(vol.est.singlePeriodFit(res2$fit, res2$DT, IntradayCol = "IntraDayBucket"))/avgXMinVol
########################### ACF  #########################
bvSDTfff[, sCorrect := simpleCorrection[IntraDayBucket],]
bvSDTfff[, LogVolsCorrect := LogVol - sCorrect,]
bvSDTfff[, nCorrect := correctionNice[IntraDayBucket],]
bvSDTfff[, LogVolnCorrect := LogVol - nCorrect,]
# bvSDTfff[, correctionOverfit := simpleCorrection[IntraDayBucket],]
# bvSDTfff[, VolOCorrect := Vol/correctionOverfit,]
par(mfrow=c(2,2))
acf(bvSDTfff$LogVol, lag.max = nIntradayBuckets*3)
acf(bvSDTfff$LogVolsCorrect, lag.max = nIntradayBuckets*3)
acf(bvSDTfff$LogVolnCorrect, lag.max = nIntradayBuckets*3)
# acf(bvSDTfff$VolOCorrect, lag.max = nIntradayBuckets*5)
par(mfrow=c(1,1))

acf(bvSDTfff$LogVol, lag.max = nIntradayBuckets)
par(mfrow=c(2,1))
acf(bvSDTfff$LogVolsCorrect, lag.max = nIntradayBuckets)
acf(bvSDTfff$LogVolnCorrect, lag.max = nIntradayBuckets)
par(mfrow=c(1,1))
plot(bvSDTfff$LogVolsCorrect[1:nIntradayBuckets*100]-bvSDTfff$LogVolnCorrect[1:nIntradayBuckets*100], type = "l")


################### Stationarity ##########################
require(fUnitRoots)
mydata <- bvSDTfff$LogVolnCorrect
mydata <- bvSDTfff$LogVolsCorrect

plot(mydata[1:100])
ppShort <- PP.test(mydata, lshort = T)
ppShort
ppLong  <- PP.test(mydata, lshort = F)
ppLong

for(i in 1:300){
adfres <- adfTest(mydata,lags = i, type = "nc")
print(c(i,adfres@test$p.value))
}
# vol.est.AICc(adfres@test$lm)
summary(adfres@test$lm)

#Always good

################ ROUGH #########################

