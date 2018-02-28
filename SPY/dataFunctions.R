if (!require(data.table)){
  install.packages("data.table")  
  require(data.table)
}

selectDays<-function(x, startDate, endDate = NA, nDays = NA){
  #Returns row of x with x$Date between startDate and endDate (or from startDate and nDays - 1 days ahead).
    #i.e. nDays = 1 yields only the start day.
  
  #x data.table with columns Date (as.date)
  #startDate as.date
  #either endDate (as.date) xor nDays (int)
  
  
  #Attempt at basic speed increase
  unDates<-unique(x$Date)
  startIndex<-match(TRUE,unDates>=startDate)
  startDate<-unDates[startIndex]

  
  print(paste0("Chosen start date: ", startDate))
  
  if(!xor(is.na(nDays),is.na(endDate))){
    stop("Specify either 'endDate' xor 'nDays'")
  }
  
  if(!is.na(nDays)){
    if(nDays <= 0 ){
      stop("nDays can't be non-positive" )
    } else {
      endDate <- unDates[startIndex + nDays - 1]
    }
  } else {
    endDate <- unDates[unDates >= endDate][1]
  }
  
  print(paste0("Chosen end date: ", endDate))
  #print(paste0("Total number of days: ")) could be calculated
    
  return(x[between(Date, startDate, endDate, incbounds = T)])
  
}


timePoints<-function(x, timeOffset = 0, lagOffset = 0, initialDelay = max(timeOffset, lagOffset), minLag = 0){
  #returns vector of indices of x corresponding to:
      #Every day starts at initialDelay (unless minLag is specified),
      #then offsets according to timeOffset or lagOffset. If using timeOffset, then rolls to nearest index
      #finally adds last point of the day
      #Note minLag is an additional condition to the first point of the day when using timeOffset. It ensures that
        #at least minLag lags are on the day before the start time is chosen.
  

  
  #x data.table with Date, Time column
  #timeOffset (double) offset in seconds 
  #lagOffset (int) offset in number of observations (referred to as lag)
  #initialDelay (double or int) how far to move into each day before choosing the first point.
    #can only specify either timeOffset or lagOffset, and assumes initialDelay 
    #to be on same "scale" i.e. either seconds or lag
  #minLag (int), explained above
  
  
  #sorts in-memory (pointer)
  setkey(x, Date)
  
  unDates<-unique(x$Date)
  dateIndicesOffset <- match(unDates, x$Date)
  if(!(timeOffset>0 | lagOffset>0)){
    stop("Specify either 'timeOffset' xor 'timeOffset'")
  }
  
  if((lagOffset>0 & minLag>0)){
    stop("minLag should only be used with timeOffset, not with lagOffset")
  }
  
  indexList<-list()
  
  #Loop across days
  for(i in 1:length(unDates)){
    currentData <- x[x$Date == unDates[i], "Time"]
    
    #timeOffset and lagOffset basically seperate functions, as lag is much simpler
    if(timeOffset > 0){
      startTime <- currentData$Time[1]
      lastTime <- currentData$Time[length(currentData$Time)]
      startTime <- startTime + initialDelay
      
      if(sum(currentData$Time<=startTime)<minLag){
        startTime <- currentData$Time[minLag]
        warning("Minimum lag not exceeded by initial delay. Initial timepoint moved according to minimum lag requirement.")
      }
      
      lengthOut<-floor((lastTime-startTime)/timeOffset)+1
      
      #use all timepoints if timeOffset is extremely small
      if(lengthOut>length(currentData$Time)){
        
        indexList[[i]]<- 1:lengthOut+dateIndicesOffset[i]-1
        warning("timeOffset so small that all timepoints are used.")
        
      } else {
      
        allTimes<-seq(startTime, by = timeOffset, length.out = lengthOut)
        
        #data.table method for finding (i.e. rolling to) indices:
        currentData[,val:=Time]
        setattr(currentData, "sorted", "Time")
        indexList[[i]]<-currentData[J(allTimes), .I, roll = "nearest", by = .EACHI][,Time:=NULL]
        
        #Add final if needed
        if(indexList[[i]][length(indexList[[i]])]!=currentData$Time[length(currentData$Time)]){
          indexList[[i]]<- rbindlist(list(indexList[[i]], list(length(currentData$Time))))+dateIndicesOffset[i]-1
        }
      }
      
    
    } else {#using lagOffset
      
      startTime <- initialDelay
      lastTime <- length(currentData$Time)
      allTimes<-seq(startTime, by = lagOffset, length.out = floor((lastTime-startTime)/lagOffset)+1)
      
      #add last, if needed
      if(allTimes[length(allTimes)]!=lastTime){
        allTimes<-c(allTimes, lastTime)
      }
      
      indexList[[i]]<-data.table(I=allTimes+dateIndicesOffset[i]-1)
      
    }
  }
  
  #returns data.table with 1 column named "I"
  return(rbindlist(indexList)$I)
  
}

prepareForEstimation<-function(x){
  #x data.table with columns: Date, Time, wprice
  
  #returns data.table:
      #log-return, convention first obs is log(wprice[2]) - log(wprice[1])
      #Time, rescaled to running from time point 0.
      #Date, same as input date. Important for ability of running timePoints()
  logRet  <- diff(log(x$wprice))
  timeCol <- x$Time[-1]-x$Time[1]
  dateCol <- x$Date[-1]
  return(data.table(Date = dateCol, Time = timeCol, Y = logRet))
  
}
