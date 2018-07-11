#### source functions neeeded
cd_UNIQUE_NAME_PER_TABLE <- getwd()

setwd(paste0(Sys.getenv("masters-thesis"), "\\Vol"))
source("vol_estimators.R")
source("BV_Analysis_Func.R")
source("FlexibleFourierForm_Func.R")
source("adf_fit_func.R")
source("PR_Func.R") # PR = Persistence and roughness


setwd(cd_UNIQUE_NAME_PER_TABLE)
rm(cd_UNIQUE_NAME_PER_TABLE)
############


dt <- BV.get_SPY_data()
bucketLengthInMinutes <- c(5,10,15,30,39,65,130,195,390)

nMinM <- nMaxM <- 100
OLS_betaArray <- array(data = NA, dim = c(nMinM, nMaxM, length(bucketLengthInMinutes)))
PVec <- rep(NA, length(bucketLengthInMinutes))
persistenceDTList <- list()
xyList <- list()
for(i in seq_along(bucketLengthInMinutes)){
  # i<-1
  print(i)
  #Compute IV given bucket lengths
  temp <- BV.data_deseason_BV_Func(dt = dt, bucketLengthInMinutes = bucketLengthInMinutes[i])
  bvSDTfff <- temp$bvSDTfff
  PVec[i] <- temp$optimalP
  
  rm(temp)

  #Quick extract/rename columns
  persistenceDT <- PR.persistence_prep_DT(bvSDTfff)
  persistenceDTList[[i]] <- persistenceDT
  
  lastBeforeNeg <- which(as.vector(acf(persistenceDT$LogVolCorrected, lag.max = length(persistenceDT$LogVolCorrected)/2, plot = F)[[1]])<1e-10)[1] - 1
  
  if(is.na(lastBeforeNeg)){
    lastBeforeNeg <- floor(length(persistenceDT$LogVolCorrected)/2)
  }

  full <- floor(seq(from = 1, to = lastBeforeNeg, length.out = nMinM+1))
  minm <- full[1:(length(full) - 1)]
  maxm <- full[2:length(full)]
  xyList[[i]] <- full
  for(j in seq_along(minm)){
    for(k in seq_along(maxm)){
      print(c(j,k))
      if(minm[j]<maxm[k]){
        OLS_betaArray[j,k,i] <- PR.est.beta(persistenceDT = persistenceDT, minm[j], maxm[k], bucketLengthInMinutes[i])
      } else {
        next
      }
    }
  }
}
# saveRDS(xyList, "xyList.rds")
# saveRDS(persistenceDTList, "persistenceDTList.rds")
# saveRDS(OLS_betaArray, "persistence_OLS_betaArray.rds")
# saveRDS(PVec, paste0("persistence_PVec__Buckets_",
                     # min(bucketLengthInMinutes), "-", max(bucketLengthInMinutes),".rds"))
# OLS_betaArray <- readRDS("persistence_OLS_betaArray.rds")
# persistenceDTList <- readRDS("persistenceDTList.rds")
# xyList <- readRDS("xyList.rds")
require(plotly)
slice <- 1
betaSlice <- OLS_betaArray[,,slice]

betaSlice2 <- betaSlice
betaSlice2[betaSlice2<quantile(betaSlice,probs = 0.05, na.rm = T)|betaSlice2>quantile(betaSlice,probs = 0.95, na.rm = T)] <- NA

rows <- xyList[[slice]][-length(xyList[[slice]])]
cols <- xyList[[slice]][-1]



p <- plot_ly(z = betaSlice, x = rows, y = cols) %>% add_surface()  %>%
  layout(
    title = "Persistence, 5-minute bucket",
    scene = list(
      xaxis = list(title = "min"),
      yaxis = list(title = "max"),
      zaxis = list(title = "&beta hat")))


p

p <- plot_ly(z = ~betaSlice2, x = rows, y = cols) %>% add_surface()
p




##############
plotLogACF <- function(persistenceDT, TradingDayLagMin, TradingDayLagMax, bucketLengthInMinutes, SPY_Bool = F){
  AnnualizedBucketLength <- bucketLengthInMinutes/(60*24*252)
  if(SPY_Bool){
    TradingDayACF <-as.vector(acf(persistenceDT$LogVolCorrected, lag.max = TradingDayLagMax, plot = F)[[1]])[TradingDayLagMin:TradingDayLagMax]
  } else {
    TradingDayACF <- PR.acfFunc(DT = persistenceDT, lags = TradingDayLagMin:TradingDayLagMax, logVolCol = "LogVolCorrected")$acf
  }
  
  if(any(TradingDayACF<0)){ # if any negative, change max, and do it again
    lastBeforeNeg <- which(TradingDayACF<0)[1] -1
    print(c(TradingDayLagMax, lastBeforeNeg))
    TradingDayLagMax <- lastBeforeNeg
    
    if(SPY_Bool){
      TradingDayACF <-as.vector(acf(persistenceDT$LogVolCorrected, lag.max = TradingDayLagMax, plot = F)[[1]])[TradingDayLagMin:TradingDayLagMax]
    } else {
      TradingDayACF <- PR.acfFunc(DT = persistenceDT, lags = TradingDayLagMin:TradingDayLagMax, logVolCol = "LogVolCorrected")$acf
    }
  }
  
  TradingDayh <- AnnualizedBucketLength*TradingDayLagMin:TradingDayLagMax
  
  plot(log(TradingDayh), log(TradingDayACF), xlab = expression(log~(h)), ylab = expression(log~ (rho)), sub = paste(bucketLengthInMinutes,"-minute buckets"))
  return(NULL)
}

plotLogACF(persistenceDTList[[slice]], xyList[[slice]][1], xyList[[slice]][length(xyList[[slice]])], bucketLengthInMinutes[slice])
PR.est.beta(persistenceDTList[[slice]], 19344^(1/4), 19344^(1/3), bucketLengthInMinutes[slice])



## BTC
persistenceDT <- PR.persistence_prep_DT(bvSDTfff = bvS_List$bvSDTfff)
plotLogACF(persistenceDT = persistenceDT, TradingDayLagMin = 1, TradingDayLagMax = length(bvSDTfff[[1]]/2), bucketLengthInMinutes = 5, SPY_Bool = T)
