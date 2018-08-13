setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/rescaling.R")
source("kernels/kernels.R")
source("spy/dataFunctions.R")
source("spy/datahandling.R")

# ESTIMATION PARAMETERS
hd <- 600 #(seconds)
hv <- hd                                      # TRY RATIO 1
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)

# GET DATA
data <-data.getFull()

# MARK DAYS
data <- data.dayID(data)

# T
Tdata <- data.TforId(data, "day", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)

# BIT PENGE

# GET DATA AND ROLL
bitfinex <- data.getbitdata("bitfinex")
bitmex <- data.getbitdata("bitmex")
kraken <- data.getbitdata("kraken")

bitfinex <- bitfinex[, all := 1]
bitmex <- bitmex[, all := 1]
kraken <- kraken[, all := 1]

bitfinex.T<-data.TforId(bitfinex, "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
bitmex.T<-data.TforId(bitmex, "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
kraken.T<-data.TforId(kraken, "all", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)

if(F){
  # TEST if works
  # fill up with 1:n per day
  #Tdata$Tval = NULL
  #ns<-Tdata[, .N, by = day]$N
  #Tval <- NULL
  #for(i in 1:length(ns)){
  #  Tval<- c(Tval,seq(1, ns[i], 1))
  #}
  #Tdata[, Tval:=Tval]
}

#################
# SETUP BUCKETS #
#################

Tdist<-function(dataT){
  # ABSOLUTELY
  absT.unsorted<-dataT$Tval
  
  # STRIP DATE
  t <- strftime(dataT$DateTime, format = "%H:%M:%S", tz = "UTC")
  
  # Set all to same date
  ntimes.unsorted <- as.POSIXct(t, format = "%H:%M:%S", tz = "UTC")
  
  # SORT
  pre.sort <-data.table(Numeric = as.numeric(ntimes.unsorted), DateTime = ntimes.unsorted, absT = absT.unsorted)
  sorted <- setkey(pre.sort)
  
  ntimes<-sorted$DateTime
  absT<-sorted$absT
  
  # regular dates
  dates<-seq(data.floor_date(ntimes[1], unit = "mins"), ntimes[length(ntimes)], by = t.freq)
  dates<-dates[-length(dates)]
  # match T date to regular date
  ind <- data.date_To_tindex(ntimes, dates)
  
  a <- data.table(DateTime = ntimes, absT = absT, index = ind)
  
  # dist
  dist <- a[, var((.SD[["absT"]])), by = index]
  
  # plot data
  return(data.table(Time = dates, Mean = dist$V1))
}

spy<-Tdist(Tdata)
finex <- Tdist(bitfinex.T)
mex <- Tdist(bitmex.T)
krak <- Tdist(kraken.T)

require(ggplot2)

buckets <- spy$Time
timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])

labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")

g1 <- ggplot() +
  geom_point(aes(x = Time, y = Mean, group = 1, col = "SPY"), data = spy) +
  ylab("Variance of the T-statistic") +
  geom_hline(yintercept = 1) +
  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
g1


#BITCOIN SECTION
if(F){
  # FINEX
  buckets <- finex$Time
  timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                  buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])
  labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")
  g2 <- ggplot() +
    geom_point(aes(x = Time, y = Mean, group = 1, col = "BitFinex"), data = finex) +
    ylab("Variance of the T-statistic") +
    geom_hline(yintercept = 1) +
    scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
  g2
  
  # MEX
  buckets <- mex$Time
  timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                  buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])
  labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")
  g3 <- ggplot() +
    geom_point(aes(x = Time, y = Mean, group = 1, col = "BitMEX"), data = mex) +
    ylab("Variance of the T-statistic") +
    geom_hline(yintercept = 1) +
    scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
  g3
  
  # KRAKEN
  buckets <- krak$Time
  timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                  buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])
  labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")
  g4 <- ggplot() +
    geom_point(aes(x = Time, y = Mean, group = 1, col = "Kraken"), data = krak) +
    ylab("Variance of the T-statistic") +
    geom_hline(yintercept = 1) +
    scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
  g4
}


########################
# TEST FOR SEASONALITY #
########################

data.Tseasonality<-function(Tdistdata, name = "Data SPY", log = F){
    # Tdistdata contains Mean | Time
    
    # Distance from mid-day
    time <- abs(as.numeric(Tdistdata$Time)-mean( c(as.numeric(Tdistdata$Time[1]), as.numeric(Tdistdata$Time[length(Tdistdata$Time)])) ))
    
    # SETUP MODEL
    test<-data.frame(Time = time/3600, Mean = Tdistdata$Mean)
    if(log){
      b<-lm(log(Mean)~Time, data = test)
    }
    else{
      b<-lm(Mean~Time, data = test)
    }
    
    print(summary(b))
    
    # MODEL PREDICTS
    if(log){
      model <- exp(as.numeric(b$coefficients[1])+as.numeric(b$coefficients[2])*test$Time)
    }
    else{
      model <- as.numeric(b$coefficients[1])+as.numeric(b$coefficients[2])*test$Time
    }
    Tdistdata$model <- model
    
    # PLOTTERIA
    buckets <- Tdistdata$Time
    timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                    buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])
    labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")
    
    g1 <- ggplot() +
      geom_point(aes(x = Time, y = Mean, col = name), data = Tdistdata) +
      geom_line(aes(x = Time, y = model, col = "Model"), data = Tdistdata, size = 1) +
      ylab("Average absolute T-statistic") +
      scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
    return(g1)
  }
  
#data.Tseasonality(mex, name = "Data bitMEX", log = F)
#data.Tseasonality(spy, name = "Data SPY", log = F)

#######################
# DIFFERENT BANDWIDTH #
#######################

# T
Tdata <- data.TforId(data, "day", hd = hd, hv = hv, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
spy<-Tdist(Tdata)

#Tdata2 <- data.TforId(data, "day", hd = hd*2, hv = hv*2, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
#spy2<-Tdist(Tdata2)

#Tdata3 <- data.TforId(data, "day", hd = hd/2, hv = hv/2, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
#spy3<-Tdist(Tdata3)

#Tdata4 <- data.TforId(data, "day", hd = hd/4, hv = hv/4, t.freq = t.freq, lag = lag, offset = offset, offset_perId = T)
#spy4<-Tdist(Tdata4)

buckets <- spy$Time
timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])

labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")

g1 <- ggplot() +
  geom_point(aes(x = Time, y = Mean, group = 1, col = "300"), data = spy) +
  #geom_point(aes(x = Time, y = Mean, group = 1, col = "600"), data = spy2) +
  #geom_point(aes(x = Time, y = Mean, group = 1, col = "150"), data = spy3) +
  #geom_point(aes(x = Time, y = Mean, group = 1, col = "75"), data = spy4) +
  labs(color = expression("Bandwidth"~mu)) +
  ylab("Variance of the T-statistic") +
  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
g1


#data.Tseasonality(spy2, name = "Data SPY", log = F)
data.Tseasonality(spy, name = "Data SPY", log = F)
#data.Tseasonality(spy3, name = "Data SPY", log = F)
#data.Tseasonality(spy4, name = "Data SPY", log = F)


############################
# SEASONALITY OF SIGMA EST #
# --------BANDWIDTHS------ #
############################

Sigmadist<-function(dataT){
  # ABSOLUTELY
  absT.unsorted<-dataT$Sigma
  
  # STRIP DATE
  t <- strftime(dataT$DateTime, format = "%H:%M:%S", tz = "UTC")
  
  # Set all to same date
  ntimes.unsorted <- as.POSIXct(t, format = "%H:%M:%S", tz = "UTC")
  
  # SORT
  pre.sort <-data.table(Numeric = as.numeric(ntimes.unsorted), DateTime = ntimes.unsorted, absT = absT.unsorted)
  sorted <- setkey(pre.sort)
  
  ntimes<-sorted$DateTime
  absT<-sorted$absT
  
  # regular dates
  dates<-seq(data.floor_date(ntimes[1], unit = "mins"), ntimes[length(ntimes)], by = t.freq)
  dates<-dates[-length(dates)]
  # match T date to regular date
  ind <- data.date_To_tindex(ntimes, dates)
  
  a <- data.table(DateTime = ntimes, absT = absT, index = ind)
  
  print(a[, .N, by = index])
  
  # dist
  dist <- a[, mean(.SD[["absT"]]), by = index]
  
  # plot data
  return(data.table(Time = dates, Mean = dist$V1))
}

# Sigma
#spy<-Sigmadist(Tdata)
#spy2<-Sigmadist(Tdata2)
#spy3<-Sigmadist(Tdata3)
#spy4<-Sigmadist(Tdata4)

#buckets <- spy$Time
#timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
#                buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])

#labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")

#g1 <- ggplot() +
#  geom_point(aes(x = Time, y = Mean, group = 1, col = "300"), data = spy) +
  #geom_point(aes(x = Time, y = Mean, group = 1, col = "600"), data = spy2) +
  #  geom_point(aes(x = Time, y = Mean, group = 1, col = "150"), data = spy3) +
  #geom_point(aes(x = Time, y = Mean, group = 1, col = "75"), data = spy4) +
#  labs(color = expression("Bandwidth"~mu)) +
#  ylab("Average absolute Sigma") +
#  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
#g1

############################
# SEASONALITY OF MU EST    #
# --------BANDWIDTHS------ #
############################

Mudist<-function(dataT){
  # ABSOLUTELY
  absT.unsorted<-dataT$Mu
  
  # STRIP DATE
  t <- strftime(dataT$DateTime, format = "%H:%M:%S", tz = "UTC")
  
  # Set all to same date
  ntimes.unsorted <- as.POSIXct(t, format = "%H:%M:%S", tz = "UTC")
  
  # SORT
  pre.sort <-data.table(Numeric = as.numeric(ntimes.unsorted), DateTime = ntimes.unsorted, absT = absT.unsorted)
  sorted <- setkey(pre.sort)
  
  ntimes<-sorted$DateTime
  absT<-sorted$absT
  
  # regular dates
  dates<-seq(data.floor_date(ntimes[1], unit = "mins"), ntimes[length(ntimes)], by = t.freq)
  dates<-dates[-length(dates)]
  # match T date to regular date
  ind <- data.date_To_tindex(ntimes, dates)
  
  a <- data.table(DateTime = ntimes, absT = absT, index = ind)
  
  print(a[, .N, by = index])
  
  # dist
  dist <- a[, var(.SD[["absT"]]), by = index]
  
  # plot data
  return(data.table(Time = dates, Mean = dist$V1))
}

# T
#spy<-Mudist(Tdata)
#spy2<-Mudist(Tdata2)
#spy3<-Mudist(Tdata3)
#spy4<-Mudist(Tdata4)

#buckets <- spy$Time
#timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
#                buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])

#labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")

#g1 <- ggplot() +
#  geom_point(aes(x = Time, y = sqrt(300)*Mean, group = 1, col = "300"), data = spy) +
#  geom_point(aes(x = Time, y = sqrt(600)*Mean, group = 1, col = "600"), data = spy2) +
#  #  geom_point(aes(x = Time, y = sqrt(150)*Mean, group = 1, col = "150"), data = spy3) +
#  geom_point(aes(x = Time, y = sqrt(75)*Mean, group = 1, col = "75"), data = spy4) +
#  labs(color = expression("Bandwidth"~mu)) +
#  ylab("Average absolute Mu") +
#  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
#g1

############################
# SEASONALITY OF BOTH EST  #
#                          #
############################
mu <- Mudist(Tdata)
sigma <- Sigmadist(Tdata)

buckets <- mu$Time
timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])

labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")

require(ggplot2)
plotdata <- data.frame(Time = mu$Time, Mu = sqrt(hd)^2*mu$Mean, Sigma = (sigma$Mean)) # mange punkter uden ændring

offset <- 0 # mean(plotdata$Sigma)-mean(plotdata$Mu)-5*10^-6
scale <-  1 #(max(plotdata$Mu)-min(plotdata$Mu))/(max(plotdata$Sigma)-min(plotdata$Sigma)) # check this later

p <- ggplot(data = plotdata, aes(x = Time)) +
  geom_point(aes(y = Mu, color = "Mu")) + #rescale this
  #geom_point(aes(y = (Sigma-offset)*scale, colour = "Sigma")) +
  geom_point(aes(y = Sigma, colour = "Sigma")) + 
  #scale_y_continuous(sec.axis = sec_axis(trans = ~./scale+offset, name = "Sigma")) +
  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels) +
  ylab("Value")

p
