 #################################################### VARIOGRAM TABLE ########################################
 PR.variogram_prep_DT <- function(bvSDTfff){
 variogramDT <- bvSDTfff[ , c("DayBucket", "IntraDayBucket", "LogVolnCorrect")]
 names(variogramDT)[!(names(variogramDT)%in%c("DayBucket", "IntraDayBucket"))] <- "LogVolCorrected"
 return(variogramDT)
 }
 
 PR.est.alpha <- function(variogramDT, m, bucketLengthInMinutes, OLS = T){
   AnnualizedBucketLength <- bucketLengthInMinutes/(60*24*252)
   TradingDayh <- AnnualizedBucketLength*(1:m)
   variogramValues <- PR.vgFunc2(variogramDT$LogVolCorrected, m)
   if(OLS){
     fit <- lm(log(variogramValues)~ 1+log(TradingDayh))
     alpha <- (fit$coefficients[2] - 1)/2  
   } else {
     startGuess <- list(NLS_c = 1, NLS_alpha = -0.3)
     fit2 <- nls(variogramValues ~ I(NLS_c*TradingDayh^(2*NLS_alpha+1)), start = startGuess, control = nls.control(maxiter = 50, tol = 1e-04, minFactor = 1/1024))
     alpha <- fit2$m$getPars()[2]
     # plot(TradingDayh, variogramValues)
   }
   
  return(alpha)
 }
 
 #################################################### PERSISTENCE TABLE ########################################
 PR.persistence_prep_DT <- PR.variogram_prep_DT
 
 
 PR.est.beta <- function(persistenceDT, TradingDayLagMin, TradingDayLagMax, bucketLengthInMinutes){
   AnnualizedBucketLength <- bucketLengthInMinutes/(60*24*252)
   TradingDayh <- AnnualizedBucketLength*TradingDayLagMin:TradingDayLagMax
   TradingDayACF <-as.vector(acf(persistenceDT$LogVolCorrected, lag.max = TradingDayLagMax, plot = F)[[1]])[TradingDayLagMin:TradingDayLagMax]
   fit <- lm(log(TradingDayACF)~ 1+log(TradingDayh))
   beta <- -fit$coefficients[2]
   
   return(beta)
 }

#################################################### BASE FUNCTIONS ########################################
################ Variogram #################
#Overly complicated method. But correct.

PR.vgFunc <- function(DT, m){
  #Returns hInd and emperical variogram values (i.e. \hat{gamma}* from bennedsen)
  totalRows <- m*length(unique(DT$DayBucket))
  resDT <- data.table(DayBucket = rep(unique(DT$DayBucket), each = m), hInd = rep(1:m, length(unique(DT$DayBucket))), count = numeric(totalRows), sumVal = numeric(totalRows))
  DT[, PR.vgDay(.SD , m, resDT, as.integer(.BY[[1]])), by = DayBucket] #Edits resDT by reference
  resFinalDT <- resDT[, list("vg"=sum(sumVal)/sum(count)), by = hInd]
  return(resFinalDT)
}

PR.vgDay <- function(DT, m, resDT, DayBucket){
  #applied to fulldata[, .SD, by = day] (i.e. should be called on every day)
    #loops across every hInd
    #Results are written by reference to resDT (should be pre-allocated)
  LogVolCorrected <- DT$LogVolCorrected
  IntraDayBucket  <- DT$IntraDayBucket
  nIntraDayBucket <- length(unique(IntraDayBucket))
  cols <- c("count", "sumVal")
  for(i in 1L:m){
    hres <- PR.vgValue(LogVolCorrected, IntraDayBucket, i, nIntraDayBucket)
    set(resDT, ((DayBucket - 1L)*as.integer(m) + i), cols, hres)
  }
}

PR.vgValue <- function(LogVolCorrected, IntraDayBucket, hInd, nIntraDayBucket){
  # One h, one day
  x <- rep(NA, nIntraDayBucket + hInd) #Creates NA values on purpose
  x[IntraDayBucket] <- LogVolCorrected
  temp   <- (LogVolCorrected-x[IntraDayBucket+hInd])^2
  sumVal <- sum(temp, na.rm = T)
  count  <- sum(!is.na(temp))
  return(list(count=count, sumVal = sumVal))
}



####### Variogram - simple method #########
PR.vgFunc2 <- function(LogVolCorrected, m){
  res <- numeric(m)
  for(i in 1:m){
    res[i] <- mean(diff(LogVolCorrected, lag = i)^2, na.rm = T) 
  }
  return(res)
}

########################## ACF ############################

PR.acfFunc <- function(DT, lags, logVolCol = "LogVolCorrected"){
  resDT <- data.table(hInd = lags, acf = numeric(length(lags)))
  LogVolCorrected <- DT[[logVolCol]]
  lagInd <- DT$lagInd
  nBucket <- max(lagInd)
  meanVal <- mean(LogVolCorrected)
  varFac <- (length(lagInd) - 1) * var(LogVolCorrected) # Variance is 1/(n-1)* sum(...). This is multiplied by (n-1) to only get the sum
  for(i in seq_along(lags)){
    print(i)
    set(x = resDT, i = i, j = "acf", value = PR.acfValue(LogVolCorrected, lags[i], nBucket, lagInd, varFac, meanVal))  
  }
  
  return(resDT)
}

PR.acfValue <- function(LogVolCorrected, hInd, nBucket, lagInd, varFac, meanVal){
  
  x <- rep(NA, nBucket + hInd) #Creates NA values on purpose
  x[lagInd] <- LogVolCorrected
  if(sum(!is.na((LogVolCorrected-meanVal)*(x[lagInd+hInd]-meanVal)))==0){
    res <- NA
  } else {
    temp1 <- (LogVolCorrected-meanVal)*(x[lagInd+hInd]-meanVal)
    temp2 <- sum(temp1, na.rm = T)
    
    # varFac2 <- (sum(!is.na(temp1))-1)*var(LogVolCorrected[!iss.na(temp1)])
    res <- temp2/varFac
  }
  
  return(res)
}




####################################### UNTESTED METHOD  ###################################

#######in any case: only valid for small lags.
# 
# PR.acfFunc <- function(DT, m){ #BASICALLY SAME AS FOR VG
#   #Returns hInd and emperical variogram values (i.e. \hat{gamma}* from bennedsen)
#   totalMean <- mean(DT$LogVolCorrected)
#   totalRows <- m*length(unique(DT$DayBucket))
#   resDT <- data.table(DayBucket = rep(unique(DT$DayBucket), each = m), hInd = rep(1:m, length(unique(DT$DayBucket))), count = numeric(totalRows), sumVal = numeric(totalRows))
#   DT[, acfDay(.SD , m, resDT, as.integer(.BY[[1]]), totalMean), by = DayBucket] #Edits resDT by reference
#   resFinalDT <- resDT[, list("vg"=sum(sumVal)/sum(count)), by = hInd]
#   return(resFinalDT)
# }
# 
# acfDay <- function(DT, m, resDT, DayBucket, totalMean){ #ALMOST SAME AS FOR VG
#   #applied to fulldata[, .SD, by = day] (i.e. should be called on every day)
#   #loops across every hInd
#   #Results are written by reference to resDT (should be pre-allocated)
#   LogVolCorrected <- DT$LogVolCorrected
#   IntraDayBucket  <- DT$IntraDayBucket
#   nIntraDayBucket <- length(unique(IntraDayBucket))
#   cols <- c("count", "sumVal")
#   for(i in 1:m){
#     hres <- PR.acfValue(LogVolCorrected, IntraDayBucket, i, nIntraDayBucket, totalMean)
#     set(resDT, ((DayBucket - 1L)*as.integer(m) + i), cols, hres)
#   }
# }
# 
# PR.acfValue <- function(LogVolCorrected, IntraDayBucket, hInd, nIntraDayBucket, totalMean){ #NOT SAME AS VG, count is the same
#   # One h, one day
#   x <- rep(NA, nIntraDayBucket + hInd) #Creates NA values on purpose
#   x[IntraDayBucket] <- LogVolCorrected
#   temp   <- (LogVolCorrected - totalMean)*(x[IntraDayBucket+hInd] - totalMean)
#   sumVal <- sum(temp, na.rm = T)
#   count  <- sum(!is.na(temp))
#   return(list(count=count, sumVal = sumVal))
# }
