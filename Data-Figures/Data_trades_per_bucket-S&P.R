# BITCOIN DATA TESTING
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

require(data.table)

# TRADES DISTRIBUTION

sp <-data.getFull()
sp <- data.trades_dist(sp, 300, "09:30:00", "16:00:00")

# S&P
plot.data <- sp
buckets <- plot.data$Time
timebreaks <- c(buckets[1], buckets[floor(length(buckets)/4)], buckets[floor(length(buckets)*2/4)],
                buckets[floor(length(buckets)*3/4)], buckets[length(buckets)])

labels <- strftime(timebreaks, format = "%H:%M", tz = "UTC")
# PLOT
g1 <- ggplot(data = plot.data) +
  geom_point(aes(x = Time, y = N, group = 1, col = "SPY")) +
  ylab("Trades per 5min bucket") +
  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
g1

####################
# TESTING FACILITY #
####################

time<-abs(as.numeric(sp$Time)-as.numeric(as.POSIXct("2018-07-04 12:45:00", tz = "UTC")))
test<-data.frame(Time = time, N = sp$N)

a<-lm(log(N) ~ Time, data = test)
summary(a)

model <- exp(as.numeric(a$coefficients[1])+as.numeric(a$coefficients[2])*test$Time)

plot.data$model <- model

g1 <- ggplot(data = plot.data) +
  geom_point(aes(x = Time, y = N, group = 1, col = "SPY")) +
  geom_line(aes(x = Time, y = model, group = 1, col = "Model")) +
  ylab("Trades per 5min bucket") +
  scale_x_datetime("Time bucket", breaks = timebreaks, labels = labels)
g1

# Without "distance to middle" trick
b <- lm(N ~ Time, data = sp)
summary(b)
