options(digits = 14, digits.secs = 5)

#### source and get data - slow
cd <- getwd()
setwd(paste0(Sys.getenv("masters-thesis"),"/SPY"))
source("dataFunctions.R")

setwd(paste0(Sys.getenv("masters-thesis-data"),"/SPY"))
fullData<-readRDS("2014_SPY_Vol_Avg.rdata")
setwd(cd)
rm(cd)
#### various data extractions
#Date format: yyyy-mm-dd

### Testing selectDays - slow

firstDay<-selectDays(fullData, as.Date("2014-01-01"), nDays = 1)
daysInRange<-selectDays(fullData, as.Date("2014-01-02"), endDate = as.Date("2014-01-03 "))
smallData<-daysInRange[1:100,]

### Testing timePoints - fast
#Used for t.index in estimation functions
plot(timePoints(daysInRange, timeOffset = 1), main="Steep slope = many observations per second")
timePoints(smallData, timeOffset = 0.1)
timePoints(smallData, lagOffset = 8, initialDelay = 68)
timePoints(smallData, lagOffset = 8, initialDelay = 67)


#### Plotting the effects - slow
plot(daysInRange$Time, daysInRange$wprice)
tIndex<-timePoints(daysInRange, timeOffset = 120, initialDelay = 600)
points(daysInRange$Time[tIndex], daysInRange$wprice[tIndex], col = "red")

#### Preparation for study - slow
readyData<-prepareForEstimation(daysInRange)
tIndex <- timePoints(readyData, timeOffset = 120, initialDelay = 600)
plot(readyData$Time, readyData$Y)
points(readyData$Time[tIndex],readyData$Time[tIndex])