library(ggplot2)
setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")


#
lag <- 10
x_es <- 0:10

y_s <- sapply(x_es, function(x) parzenkernfunction(x,lag))

x_axis <- x_es/lag

my_df <- data.frame(cbind(x_axis,y_s))

ggplot(my_df, aes(x_axis,y_s)) + 
  geom_point(color = "blue") + 
  geom_line(color = "black") +
  xlab("") + ylab("") +
  ggtitle(TeX("The parzen kernel $w(x)$")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  
  
