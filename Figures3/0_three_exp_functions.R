library(ggplot2)
library(latex2exp)
library(reshape2) #Til Melt function


mat <- 6.5
h_n <- c(1,1/2,1/4,1/12) #1hour, 30min, 15min, 5min
time_point <- 4 #4 hours into the day

#Create X'es
all_timepoints <- seq(from = 0,to = 6.5, length.out = 23400)

#Weight function
weight_func <- function(h_n) {
  weights <- numeric(length(all_timepoints))
  weights[all_timepoints > time_point] <- 0
  weights[all_timepoints <= time_point] <- 1/h_n*exp((all_timepoints[all_timepoints <= time_point]-time_point)/h_n)
  return(weights)
}

#Create weights
weights <- sapply(h_n, function(h) weight_func(h))
colnames(weights) <- c("1hour", "30min", "15min", "5min")  

#Create data_frame for plots - with melt function or similar to understand it
plot_frame <- melt(data = weights, 
                   measure.vars = c("1hour", "30min", "15min", "5min"))
plot_frame <- cbind(all_timepoints,plot_frame[,c(2,3)])
colnames(plot_frame) <- c("time_point", "Bandwidth", "weight")

#Plot
ggplot(plot_frame, aes(time_point,weight,color = Bandwidth)) +
  geom_line() +
  xlab("Time in hours") +
  ylab(TeX('Value of $\\frac{1}{h_n} \\cdot K\\left( \\frac{t_i-t}{h_n} \\right)'))
  #ggtitle(TeX('Weights of $\\Delta X_{t_i}$ for different bandwidths at t_i = 4 hours')) +
  #theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  
  
  
