# DATA HANDLING OPTIMIZATION

setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
source("spy/dataFunctions.R")
source("spy/datahandling.R")

# LET US TAKE A LOOK AT SECTION 5.3 - REVERSION AFTER DRIFT BURSTS

# The idea is that we can 'model' the return up to X-min after the burst peak
# by looking at the return X-min down to Y-min before
# Ideally, we could expand this into something more elaborate and ~cool~
# We CANNOT use T in the model unless we calculate T way more often... (Once per second)
# OR we could calculate Ts whenever we need them - on the fly inside function...?

# FETCH DATA AND LABEL DAYS
# by labeling days, we 
data<-data.getFull()
data<-data.dayID(data, id = "day")

require(Rcpp)
require(microbenchmark)

dt <- data[day <= 2]




data.time_To_tindex<-function(times_seq, t_times){
  # ASSUME BOTH ARE ORDERED
  a <- sapply(times_seq, function(x) which.max(t_times[t_times<=x]))
  return(do.call("rbind", a))
}

data.date_To_tindex_dt<-function(dates, t_dates){
  # ASSUMES SORTED
  date <- t_dates
  sorter <- data.table(date, val = date)
  setattr(sorter, "sorted", "date")
  table <- sorter[J(dates), roll = "nearest"]
  out <- match( table$val, t_dates)
  return(out)
}

data.date_To_tindex<-function(dates, t_dates){
  # dates are the dates you wish to find the t.index for
  # t_dates are the dates that the t.index will match
  # ASSUMES ORDERED
  a <- sapply(dates, function(x) which.max(t_dates[t_dates<=x]))
  return(do.call("rbind", a))
}

Date <- dt$DateTime
Time <- dt$Time
t.freq = 5


# SETUP T-INDEX #
dates<-seq(data.floor_date(Date[1]), Date[length(Date)], by = t.freq)
times_seq <- seq(floor(Time[1]), Time[length(Time)], by = t.freq )

times <- data.date_To_tindex(dates = dates, t_dates = Date)

times2 <- data.time_To_tindex(times_seq, Time)

times3 <- data.date_To_tindex_dt(dates, Date)

dt[times,]
dt[times3,]

length(times3)-length(unique(times3))


microbenchmark(data.date_To_tindex(dates = dates, t_dates = Date),
               data.time_To_tindex(times_seq, Time),
               data.date_To_tindex_dt(dates, Date),
               times = 1)


