################ Additional data.table functions (DEALING WITH BUCKETS (for pre-averging))
vol.est.addIntradayBucket <- function(DT, nDailyBuckets){
  #Requires data.table with col "id"
  #Adds (by reference) new col to DT
  DT[, "IntradayID" := id %% nDailyBuckets]
  DT[IntradayID == 0, "IntradayID" := nDailyBuckets]
}

vol.est.SimIntradayBucket <- function(m, N, t, ...){
  bucketWidth = m * 60
  bucketIndices = seq(bucketWidth, N, by = bucketWidth) + 1
  buckets <- c(-1, t[bucketIndices]) #NOTE THE -1 is IMPORTANT (need a value < 0, to make binning work properly)
  return(list(bucketWidth = bucketWidth, buckets = buckets, bucketIndices = bucketIndices))
}

vol.est.DataIntradayBucket <- function(DT, m){
  #DT data.table with col "Date"
  if(390 %% m != 0){
    stop("m needs to specify equal-sized buckets from 9.30 to 16.00, i.e.")
  }
  
  allDays <- as.POSIXct.Date(as.Date(unique(DT$Date))) 
  
  buckets2 <- seq(from = 34200+m*60, to = 57600, by = m*60) #every day, get hours 09.30+m to 16.00
  buckets <- .POSIXct(character(length(allDays)*length(buckets2)), tz = "UTC") #rep(NA, length(allDays)*length(buckets2))
  
  for(i in seq_along(allDays)){
    buckets[(((i-1)*length(buckets2)):(i*length(buckets2)))[-1]] <- allDays[i]+buckets2+0.0001
  }
  
  return(buckets)
}
  
################ IV estimation - data.table functions
require(data.table)

##### Reduce from mutiple paths to final estimates
vol.est.IVest <- function(DT, cols = "logPrice", maxId){
  #Wrapper for vol.est.IVestPathwise.
  
  # INIT
  M <- length(cols)
  if(missing(maxId)){
    maxId <- last(DT$id)
  }
  
  # ESTIMATES
  PathwiseRes <- vol.est.IVestPathwise(DT, cols, maxId)
  
  # COLLAPSE PATHWISE RESULTS
  meansDT <- data.table(id = 1:maxId, matrix(rep(NA, maxId * 4), ncol = 4))
  names(meansDT)[-1] <- c("rv", "bv", "rvS", "bvS")
  
  if(M > 1){
    for(i in seq_along(PathwiseRes$IVest)){
      meansDT[,i+1] <- rowMeans(PathwiseRes$IVest[[i]][,-1])
    }
  } else {
    for(i in seq_along(PathwiseRes$IVest)){
      meansDT[, i+1] <- PathwiseRes$IVest[[i]][, -1]
    }
  }
  
  
  return(meansDT)
}

##### Function combining all below!
vol.est.IVestPathwise <- function(DT, cols = "logPrice", maxId){
  #data.table with cols
    #id
    #logPrice
  #Returns omega2 and IV estimates within rows with same id
  
  ##INIT
  #Number of paths
  M <- length(cols)
  
  if(missing(maxId)){
      maxId <- last(DT$id)
  }
  
  ##simple estimators
  rDT  <- vol.est.rByID(DT, cols)
  simpleIVList <- vol.est.simpleIV(rDT, M)
  
  ## adv estimators
  #Pre-avg
  temp <- vol.est.preAByID(DT, cols)
  rSDT <- temp$rSDT
  Kvec <- temp$Kvec
  
  #Clear memory
  rm(temp)
  gc()
  
  #omega2
  omega2est <- vol.est.omega2ByID(rDT, M)
  
  #estimates
  advIVList <- vol.est.advIV(rSDT, omega2est, Kvec, M)
  
  countDT <- DT[, .N, by = id]
  
  #List of 2, (IVest is a list of 4)
  return(list(IVest = c(simpleIVList, advIVList), omega2est = omega2est, countDT = countDT))
}

vol.est.simpleIV <- function(rDT, M = 1){
  #data.table with cols
    #id
    #rx of log-returns

  rvDT <- rDT[, lapply(.SD, function(x){sum(x^2)}), by = id]
  bvDT <- rDT[, lapply(.SD, function(x){vol.est.BV(x)}), by = id]
  
  names(rvDT)[-1] <- paste0("rv", 1:M)
  names(bvDT)[-1] <- paste0("bv", 1:M)
  
  return(list(rvDT=rvDT, bvDT=bvDT))
}

vol.est.omega2ByID <- function(rDT, M = 1){
  #data.table with cols
    #id
    #rx of log-returns
  #Vector omega2 estimates - NOTE: "id" col is removed on purpose
  return(vol.est.AddmissingIDs(rDT[, lapply(.SD, function(x){vol.est.omega2(x)}), by = id]))
}

vol.est.advIV <- function(rSDT, omega2estDT, Kvec, M = 1){
  #data.table with cols
    #id
    #rS of pre-averaged log-prices (i.e. pre-avg versions of log-returns)
  #vector of omega2 estimates w/ length = No unique id
  #vector of K's used for preaveraging w/ length = No unique id
  
  # rvSDT <- rSDT[, list(rvS = vol.est.RVstar(.SD$rS, K = Kvec[.BY$id], omega2Est = omega2estDT[.BY$id])), by = id]
  # bvSDT <- rSDT[, list(bvS = vol.est.BVstar(.SD$rS, K = Kvec[.BY$id], omega2Est = omega2estDT[.BY$id])), by = id]
  temp <- rSDT[,Kvec[.BY$id], by = id]
  rvSDT <- rSDT[, mapply(vol.est.RVstar,  paData_subset = .SD, omega2Est = omega2estDT[.BY$id], K = Kvec[.BY$id], SIMPLIFY = F), by = id]
  bvSDT <- rSDT[, mapply(vol.est.BVstar,  paData_subset = .SD, omega2Est = omega2estDT[.BY$id], K = Kvec[.BY$id], SIMPLIFY = F), by = id]
    
  names(rvSDT)[-1] <- paste0("rvS", 1:M)
  names(bvSDT)[-1] <- paste0("bvS", 1:M)
  
  return(list(rvSDT = rvSDT, bvSDT = bvSDT))
}

vol.est.preAByID <- function(DT, cols = "logPrice"){
  #Takes data.table with columns:
    #id
    #logPrices
  
  # tempDT <- DT[, 2*floor(sqrt(.N)/2), by = id]
  # #Vector with length = No. unique ID w/ data, and the ids
  # Kvec <- tempDT$V1
  # idVec <- tempDT$id
  # rm(tempDT)
  # gc()
  # 
  # missingID <- which( !(1:max(idVec) %in% idVec)) #Finds which id's at not present
  # for(i in seq_along(missingID)){ #Slow but fastest alternative so far
  #   Kvec <- append(Kvec, 0, after = missingID[i]-1)
  # }
  # # DT[, .N, by = id]
  Kvec <- vol.est.AddmissingIDs(DT[, 2*floor(sqrt(.N)/2), by = id])
  #data.table with cols: id, rSx
  # DT[, print(.BY$id), by = id]
  rSDT <- DT[, lapply(.SD, function(x) {vol.est.preA(x,  Kvec[.BY$id])}), by = id, .SDcols = cols]  #Utilizes Kvec having "full length" even if some buckets are empty
  names(rSDT)[-1] <- paste0("rS", 1:length(cols))
  
  return(list(Kvec = Kvec, rSDT = rSDT))
}

vol.est.AddmissingIDs<-function(DT){ # SLOW - but should only be used on "small" data.tables. 
  #DT data.table with columns:
    #id - has to be the first column
    #Any ONE-other column
  
  vec <- DT[[2]] # Fast extract of column 2 (non-id)
  idVec <- DT[[1]] # Fast extract of column 2 (id)

  missingID <- which( !(1:max(idVec) %in% idVec)) #Finds which id's at not present
  for(i in seq_along(missingID)){ 
    vec <- append(vec, NA, after = missingID[i]-1)
  }
  return(vec)
}

vol.est.rByID <- function(DT, cols = "logPrice"){
  #Takes data.table with columns:
    #id
    #logPrices
  
    retDT <- DT[, lapply(.SD, diff), by = id, .SDcols = cols]
    names(retDT)[-1] <- paste0("r", 1:length(cols))
  return(retDT)
}

####### STATS

vol.est.rowVar <- function(DT, rowMean){
  #Requires data.table:
    # Without id
    # With strictly more than 1 additional column
   
  #NOTE: -2 = -1 + - 1:     -1 col as we remove "id", -1 col for unbias, 
  return(rowSums((DT- rowMean)^2)/(dim(x)[2] - 1) )
}
########## BASE FUNCTIONS

vol.est.RVstar<-function(paData_subset, K, omega2Est, theta=1){
  N <- length(paData_subset) + K - 1
  phiK <- vol.est.phiK(K)
  # print(c((N/(N-K+2)),(1/(K*phiK)) , sum(paData_subset^2) ,omega2Est/(theta*phiK)))
  return((N/(N-K+2)) * (1/(K*phiK)) * sum(paData_subset^2) - omega2Est/(theta^2*phiK))
}
vol.est.BV<-function(paData_subset){
  lenD <- length(paData_subset)
  (lenD/(lenD-1)) * (pi/2) * sum(abs(paData_subset[1:(lenD-1)]*paData_subset[2:lenD]))
}

vol.est.BVstar<-function(paData_subset, K, omega2Est, theta=1){
  lenD <- length(paData_subset)
  N <- lenD + K - 1
  phiK<-vol.est.phiK(K)
  return((N/(N-2*K+2)) * (1/(K*phiK)) * (pi/2) 
         * sum(abs(paData_subset[1:(lenD-K)]*paData_subset[(K+1):lenD])) 
         - omega2Est/(theta^2*phiK))
}

vol.est.omega2 <- function(Data_subset){
  lenD <- length(Data_subset)
  return(-mean(Data_subset[1:(lenD-1)]*Data_subset[2:lenD]))
}

vol.est.phiK<-function(K){
  return((1+2*(K^(-2)))/12)
}

vol.est.preA<-function(logPrice_Data, K){
  if(!((K%%2==0)&(K>1))){
    stop("Needs K even")
  }
  
  N <- length(logPrice_Data) 
  lOut <- N-K+1
  r <- rep(NA, lOut)
  
  for(i in 1:(lOut)){
    r[i] <- sum(logPrice_Data[i+(K/2):(K-1)])-sum(logPrice_Data[i+(0):(K/2-1)])
  }  
  
  if(any(is.na(r))){
    print(1)
  }
  return(r/K)
}


