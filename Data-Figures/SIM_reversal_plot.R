setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("kernels/kernels.R")
source("estimation/estimates_reloaded.R")
source("estimation/estimates_revolution.R")

# CHECK FOR SEASONALITY OF T IN UN-EVEN DT STUDY TO DETERMINE ORIGIN OF SEASONALITY

# SIMULATE
setting <- sim.setup(Npath = 1000, Nsteps = 23400, omega = 1.6*10^-5)

set.seed(2342)

hest <- sim.heston(setting)

alpha <- 0.65
beta <- 0.1
c_1 <- (1-alpha)*0.005/(10/(60*24*7*52))^(1-alpha)
c_2 <- sqrt((1-2*beta)*(0.00093*0.25)^2/(10/(60*24*7*52))^(1-2*beta))

hest <- sim.addvb(hest, burst_time = 0.5, interval_length = 0.05, c_2 = c_2, beta = beta, reverse = F)
hest <- sim.adddb(hest, burst_time = 0.5, interval_length = 0.05, c_1 = c_1, alpha=alpha, reverse = F)

heston <- hest

hest$Y <- t(diff(t(as.matrix(hest$Y))))

# ESTIMATION PARAMETERS
hd <- 300/(52*7*24*60*60) #(seconds)
hv <- 12*hd                                          
lag = 10
t.freq = 5 # every 5 seconds
offset = 12 # skips the first minute (5*12seconds)
#t.index <- seq(offset, setting$Nsteps, by = t.freq)

mufull <- est.mu.mat.2.0(data = hest, hd = hd, bandwidth_rescale = T)#,t.index = t.index)$mu
mu <- sqrt(hd)*mufull$mu#[,]
sigma2 <- est.sigma.mat.2.0(data = hest, hv = hv, lag = lag, bandwidth_rescale = T)$sig#[,]#,t.index = t.index,lag = lag)$sig
Tstat <- mu/sqrt(sigma2)


# SINGLE PATH PLOT
mid <- 23400/2
high <- mid+120
low <- mid-120

plot(heston$Y[1,10000:15000], type = "l")
plot(Tstat[1,10000:15000], type = "l")

plot(heston$Y[1,low:high], type = "l")
plot(Tstat[1,low:high], type = "l")

#burststart <- 0.45*23400-low
burstend <- 0.50*23400-low

index <- 2
burstpeak <- which.max(abs(Tstat[index,low:high]))


times <- as.POSIXct("2018-08-03 09:30:00", "UTC")

require(ggplot2)
plotdata <- data.frame(Time = times+heston$time[low:high]*60*60*24*7*52, Price = exp(heston$Y[index,low:high]), T.statistic = Tstat[index,low:high])

offset <- mean(plotdata$Price)
scale <- (max(plotdata$T.statistic)-min(plotdata$T.statistic))/(max(plotdata$Price)-min(plotdata$Price)) # check this later

p <- ggplot(data = plotdata, aes(x = Time)) +
  geom_area(aes(y = T.statistic)) + #rescale this
  geom_line(aes(y = (Price-offset)*scale, colour = "Price"), size = 1) + # size might change
  scale_y_continuous(sec.axis = sec_axis(trans = ~./scale+offset, name = "Price")) +
  #geom_vline(xintercept = plotdata$Date[burststart], alpha = 0.4) + 
  geom_vline(xintercept = plotdata$Time[burstend], alpha = 1) +
  geom_vline(xintercept = plotdata$Time[burstpeak], alpha = 1, colour = "dodgerblue3")
p

# MINI STUDY

Rp <- heston$Y[,(11700+300)] - heston$Y[, 11700]
Rm <- heston$Y[, 11700] - heston$Y[,(11700-300)]

lm(Rp~Rm)

ind <- 1
peak <- which.max(abs(Tstat[ind,10000:20000]))
Rp <- heston$Y[ind,(peak+300)] - heston$Y[ind, peak]
Rm <- heston$Y[ind, peak] - heston$Y[ind,(peak-300)]


