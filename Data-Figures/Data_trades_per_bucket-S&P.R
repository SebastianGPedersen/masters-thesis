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

sp<-data.getFull()
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
