#Parameter. Corresponding to original T-estimator
c <- sqrt(2)
q <- qnorm(0.975)
sigma <- sqrt(0.0457)

#test
k <- seq(0,2, by = 0.001)

#Rejection percentage above
n_sims <- 10000
rejections <- numeric(length = length(k))

for (i in 1:n_sims){
  norms <- rnorm(length(k),mean = 0,sd = sigma)
  T_s <- c*k/sqrt(k^2+sigma^2) + 1/sqrt(k^2+sigma^2)*norms
  rejections <- rejections + (T_s > q)
}

rejection_perc <- rejections / n_sims
plot(k,rejection_perc*100)


#Theoretical max
theoretical_max <- 1-pnorm(sqrt(q^2-c^2))
theoretical_max*100



## My calculations were correct - the theoretical maximum is the actual maximum. /Sebastian

