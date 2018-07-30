library(ggplot2)
library(latex2exp)
library(dplyr)

k_n <- 8 #Antal med positiv vægt

kernel_function <- function(k_n,begin_index,all_time_points) {
  x_es <- seq(1:(k_n-1)) / k_n
  g_s <- sapply(x_es, function(x) min(x,1-x))
  
  weights <- numeric(length(all_time_points))
  weights[begin_index+1:length(x_es)] <- g_s
  return(weights)
}

### PLOTS

### All the individual
time_points <- 0:30
begin_points <- 1:(length(time_points)-k_n)
first_try <- kernel_function(k_n, 1, time_points)

tmp <- list()
for (i in 1:length(begin_points)) {
  tmp[[i]] <- data.frame(
    x_akse = time_points,
    y_akse = kernel_function(k_n, begin_points[i], time_points),
    grp = i,
    trans = 0.3+0.7*i/length(begin_points)
  )
}

plot_frame <- do.call(rbind, tmp)


### The full - group by sum
all_points <- plot_frame %>% group_by(x_akse) %>% summarise(full_sum = sum(y_akse))

#ggplot
ggplot() +
  geom_line(data = plot_frame, aes(x_akse,y_akse, group = grp, alpha = I(trans), color = "a")) +
  geom_line(data = all_points, aes(x_akse, full_sum, color = "b")) +
  xlab("Time points") + ylab("Weight") +
  ggtitle(TeX('Weight of $\\Delta Y_i$')) +
  theme(plot.title = element_text(hjust = 0.5, size = 14)) +
  scale_color_manual(name = "Functions", values = c("black", "blue"),
                     labels = unname(TeX(
                       c("Pre-avg. functions",' Sum of weights'))))


