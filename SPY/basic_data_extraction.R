options(digits = 14, digits.secs = 5)

if (!require(data.table)){
  install.packages("data.table")  
  require(data.table)
}

setwd(paste0(Sys.getenv("master-thesis"),"/SPY"))
getwd()

source("dataFunctions.R")


# setwd(.../Dropbox/Speciale/Data/SPY)
test<-readRDS("2014_SPY_Vol_Avg.rdata")



test2<-selectDays(test, as.Date("2014-01-01"), nDays = 1)
test4<-selectDays(test, as.Date("2014-01-02"), endDate = as.Date("2014-01-03 "))
test3<-test2[1:100,]

plot(timePoints(test2, timeOffset = 1)$I, main="Steep slope = many observations per second")
sum(diff(timePoints(test4, timeOffset = 1)$I)==0)
timePoints(test3, timeOffset = 0.1)
timePoints(test3, lagOffset = 8, initialDelay = 68)
timePoints(test3, lagOffset = 8, initialDelay = 67)



