options(digits = 14, digits.secs = 5)

#### source functions neeeded
cd_UNIQUE_NAME_SPY_IV_ <- getwd()

setwd(paste0(Sys.getenv("masters-thesis"), "\\Vol"))
source("vol_estimators.R")
source("FlexibleFourierForm_Func.R")
source("adf_fit_func.R")

setwd(cd_UNIQUE_NAME_SPY_IV_)
rm(cd_UNIQUE_NAME_SPY_IV_)


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
dt<-dt[dt$Date != "2014-10-30"]

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
resListIntraday <- vol.est.IVestPathwise(dt)
# resListIntraday <- readRDS("5minresListIntraday.rds")
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
resListIntraday
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
res <- FFF.estOptimalP(DT = bvSDTfff, maxP =  nIntradayBuckets/2 - 1, dataCol = "LogVol", IntradayCol = "IntraDayBucket", dayCol = "DayBucket", dailyInteraction = F, dailyParameter = F)
res$CritVec
res$fit$rank
summary(res$fit)
sum(is.na(res$fit$coefficients))
res$optimalP
res$fit

plot(res$DT$LogVol[1:nIntradayBuckets])
lines(res$fit$fitted.values[1:nIntradayBuckets])
plot(res$fit$fitted.values[1:nIntradayBuckets]-res$fit$coefficients[1])

plot(exp(FFF.singlePeriodFit(res$fit, res$DT, IntradayCol = "IntraDayBucket")), type = "l")

##### Alternate simple way - better for Jensen's inequality 
### Taking log first to bring in line with fff
fit <- lm(bvSDTfff$LogVol ~ 1 + as.factor(bvSDTfff$IntraDayBucket))
plot(FFF.singlePeriodFit(fit = fit,DT = bvSDTfff, IntradayCol = "IntraDayBucket"), type = "l", col = "blue")
simpleCorrection <- FFF.singlePeriodFit(fit = fit,DT = bvSDTfff, IntradayCol = "IntraDayBucket")

##### PLOT
plot(avgIntradayVol$IntraDayBucket, simpleCorrection, type = "l", lwd = 2, xlab = "s", ylab = "Estimated log volatility", 
                            main = "Estimated intraday seasonality", sub = paste0(bucketLengthInMinutes, "-minute buckets. P = ", res$optimalP))
correctionNice <- FFF.singlePeriodFit(res$fit, res$DT, IntradayCol = "IntraDayBucket")
lines(correctionNice, type = "l", col = "red", lwd = 2)
legend("topright", legend=c("Simple", "Flexible Fourier form"),
       col=c("black",  "red"), lty=1, cex=0.8)
# lines(exp(FFF.singlePeriodFit(res2$fit, res2$DT, IntradayCol = "IntraDayBucket"))/avgXMinVol, type = "l", col = "blue")
# correctionOverfit <- exp(FFF.singlePeriodFit(res2$fit, res2$DT, IntradayCol = "IntraDayBucket"))/avgXMinVol
########################### ACF  #########################
bvSDTfff[, sCorrect := simpleCorrection[IntraDayBucket],]
bvSDTfff[, LogVolsCorrect := LogVol - sCorrect,]
bvSDTfff[, nCorrect := correctionNice[IntraDayBucket],]
bvSDTfff[, LogVolnCorrect := LogVol - nCorrect,]
IntradayPerDay <- (60/bucketLengthInMinutes)*24
bvSDTfff[, lagInd:= (as.numeric(DayBucket)-1)*IntradayPerDay + (IntraDayBucket)]
# bvSDTfff[, correctionOverfit := simpleCorrection[IntraDayBucket],]
# bvSDTfff[, VolOCorrect := Vol/correctionOverfit,]

## Bitcoin
bvSDTfff <- bvS_List$bvSDTfff
nIntradayBuckets <- max(bvSDTfff$IntraDayBucket)
IntradayPerDay <- (60/bucketLengthInMinutes)*24
##

lags <- 1:(nrow(bvSDTfff)/2)
acfNoCorrectionRes <- PR.acfFunc(bvSDTfff, lags, logVolCol = "LogVol")
acfRes <- PR.acfFunc(bvSDTfff, lags, logVolCol = "LogVolnCorrect")
acfRes2 <- PR.acfFunc(bvSDTfff, lags, logVolCol = "LogVolsCorrect")
# temp <- bvSDTfff[1:78]
# PR.acfFunc(temp, 1:5, logVolCol = "LogVol")
# acf(x = temp$LogVol, lag.max = 5, plot = F)$acf

# plot(acfNoCorrectionRes[1:(IntradayPerDay*3)], ylab = "ACF", xlab = "Lag in buckets", main = "ACF of Log Volatility - Alternate modelling", sub = "17.5 hour night")
# barplot(height = acfNoCorrectionRes[1:(IntradayPerDay*3)]$acf , names.arg = acfNoCorrectionRes[1:(IntradayPerDay*3)]$hInd, space = 0, ylab = "ACF", xlab = "Lag in buckets", main = "ACF of Log Volatility - Alternate modelling", sub = "17.5 hour night")
# barplot(height = acfRes[1:(IntradayPerDay*3)]$acf , names.arg = acfRes[1:(IntradayPerDay*3)]$hInd, space = 0, ylab = "ACF", xlab = "Lag in buckets", main = "ACF of deseasonalized Log Volatility - Alternate modelling", sub = "17.5 hour night")

### Bitcoin
## Effect of deseasonalizing - none
par(mfrow = c(3,1))
maxLagFact <- min(IntradayPerDay*3, max(lags))
barplot(height = acfNoCorrectionRes[1:maxLagFact]$acf , names.arg = acfNoCorrectionRes[1:maxLagFact]$hInd, space = 0, ylab = "ACF", xlab = "Lag in buckets", main = "ACF of Log Volatility - Alternate modelling", sub = "17.5 hour night")
barplot(height = acfRes[1:maxLagFact]$acf , names.arg = acfRes[1:maxLagFact]$hInd, space = 0, ylab = "ACF", xlab = "Lag in buckets", main = "ACF of deseasonalized Log Volatility - Alternate modelling", sub = "17.5 hour night")
barplot(height = acfRes2[1:maxLagFact]$acf , names.arg = acfRes2[1:maxLagFact]$hInd, space = 0, ylab = "ACF", xlab = "Lag in buckets", main = "ACF of deseasonalized Log Volatility - Alternate modelling", sub = "17.5 hour night")
par(mfrow = c(1,1))
any(diff(acfRes[1:maxLagFact]$hInd) != 1)

## long acf
par(mfrow = c(3,1))
maxLagFact <- min(5000, max(lags))
barplot(height = acfNoCorrectionRes[1:maxLagFact]$acf , names.arg = acfNoCorrectionRes[1:maxLagFact]$hInd, space = 0, ,ylab = "ACF", xlab = "Lag", main = "ACF of Log Volatility")
barplot(height = acfRes[1:maxLagFact]$acf , names.arg = acfRes[1:maxLagFact]$hInd, space = 0, ylab = "ACF", xlab = "Lag", main = "ACF of deseasonalized Log Volatility")
barplot(height = acfRes2[1:maxLagFact]$acf , names.arg = acfRes2[1:maxLagFact]$hInd, space = 0, ylab = "ACF", xlab = "Lag", main = "ACF of deseasonalized Log Volatility")
par(mfrow = c(1,1))
any(diff(acfRes[1:maxLagFact]$hInd) != 1)
###

require(ggplot2)
q <- ggplot(data = acfNoCorrectionRes[1:(IntradayPerDay*3)], mapping = aes(x = hInd, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = hInd, yend = 0))+ 
  
  theme(panel.grid.major = element_blank(), 
                                                            panel.grid.minor = element_blank(),
                                                            panel.background = element_blank(), 
                                                            axis.line = element_line(colour = "black"),  
                                                            panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  coord_cartesian(ylim = c(0.05, 1)) +
  ggtitle("Alternate modelling") +  ylab("ACF") + xlab("Lag in buckets") +
  theme(plot.title = element_text(hjust = 0.5))





vals <- acf(bvSDTfff$LogVol, lag.max = nIntradayBuckets*3, xlab = "Lag in buckets", main = "ACF of Log Volatility - Chosen modelling",  sub = "0 hour night", ci = 0)
regularDF <- data.frame(hInd = vals$lag, acf = vals$acf)[-1,]
q2 <- ggplot(data = regularDF, mapping = aes(x = hInd, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = hInd, yend = 0))+ 
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),  
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  coord_cartesian(ylim = c(0.05, 1)) +
  ggtitle("Chosen modelling") +  ylab("ACF") + xlab("Lag in buckets") +
  theme(plot.title = element_text(hjust = 0.5))
q2

# install.packages("gridExtra")
library(gridExtra)
grid.arrange(q2, q, ncol=2, nrow = 1)

par(mfrow=c(2,1))
plot(acfNoCorrectionRes[1:(IntradayPerDay)], xlab = "lags", main = "ACF for log-vol")
plot(acfRes[1:(IntradayPerDay)])
par(mfrow=c(1,1))

acf(bvSDTfff$LogVol, lag.max = length(bvSDTfff$LogVol)/2)
acf(bvSDTfff$LogVolnCorrect, lag.max = length(bvSDTfff$LogVolnCorrect)/2)

par(mfrow=c(3,1))
acf(bvSDTfff$LogVol, lag.max = nIntradayBuckets*3, main = "Log volatility", ci = 0)
acf(bvSDTfff$LogVolsCorrect, lag.max = nIntradayBuckets*3, main = "Simple deseasonalized log volatility", ci = 0)
acf(bvSDTfff$LogVolnCorrect, lag.max = nIntradayBuckets*3, main = "Flexible Fourier form deseasonalized log volatility", ci = 0)
# acf(bvSDTfff$VolOCorrect, lag.max = nIntradayBuckets*5)
par(mfrow=c(1,1))

acf(bvSDTfff$LogVol, lag.max = nIntradayBuckets)
par(mfrow=c(2,1))
acf(bvSDTfff$LogVolsCorrect, lag.max = nIntradayBuckets)
acf(bvSDTfff$LogVolnCorrect, lag.max = nIntradayBuckets)
par(mfrow=c(1,1))
plot(bvSDTfff$LogVolsCorrect[1:nIntradayBuckets*100]-bvSDTfff$LogVolnCorrect[1:nIntradayBuckets*100], type = "l")

# saveRDS(bvSDTfff, "asd")
################### Stationarity ##########################
require(fUnitRoots)
bvSDTfff <- bvS_List$bvSDTfff
mydata <- bvSDTfff$LogVolsCorrect
mydata <- bvSDTfff$LogVolnCorrect


plot(mydata[1:100])
ppShort <- PP.test(mydata, lshort = T)
ppShort
ppLong  <- PP.test(mydata, lshort = F)
ppLong
maxlag <- 1
maxADF <- numeric(maxlag)

for(i in 1:maxlag){
  maxADF[i] <- adfTest(mydata,lags = i, type = "nc")@test$p.value
  print(c(i, maxADF[i]))
}

# adfres <- adfTest(mydata,lags = 20, type = "nc")
# FFF.AICc(adfres@test$lm)
summary(adfres@test$lm)

max(maxADF)
adfTest(mydata,lags = 65, type = "nc")@test$p.value
#Always good
pacf(mydata, lag.max = 1000)



### BC
bvSDTfff <- bvS_List$bvSDTfff
mydata <- rep(NA, max(bvSDTfff$lagInd))
mydata[bvSDTfff$lagInd] <- bvSDTfff$LogVolnCorrect
### 

resList <- adf.Fit.optimalP(y = mydata, maxP = 300, critF = FFF.AICc)
resList$optimalP
# summary(resList$fit)
#### SPY RESULTS
# 5-min: 162
# 10-min: 37
# 15-min: 52
# 30-min: 16
plot(resList$CritVec)
resList$NcompleteCases
plot(resList$testStats)
max(resList$testStats[resList$NcompleteCases>=300 & !is.na(resList$testStats)])
plot(resList$testStats[!is.na(resList$testStats)])
plot(resList$testStats[resList$NcompleteCases>=400])
# resList30 <- resList
# resList <-  resList5
# resList <-  resList10



plot(1:max(bvSDTfff$IntraDayBucket), bvSDTfff$sCorrect[1:max(bvSDTfff$IntraDayBucket)], type = "l", lwd = 2, xlab = "s", ylab = "Estimated log volatility", 
     main = "Estimated intraday seasonality", sub = paste0(bucketLengthInMinutes, "-minute buckets. P = ", res$optimalP))

lines(bvSDTfff$nCorrect[1:max(bvSDTfff$IntraDayBucket)], type = "l", col = "red", lwd = 2)
legend("topright", legend=c("Simple", "Flexible Fourier form"),
       col=c("black",  "red"), lty=1, cex=0.8)


plot()
lines(, col = "red")