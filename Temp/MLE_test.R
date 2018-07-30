X_es <- rnorm(n =1000, mean = 0,sd = 0.01)
x_t <- X_es[2:1000]
x_t_1 <- X_es[1:999]
n <- length(x_t)

rhos <- seq(-0.9,0.9, by = 0.005)

#Ren log function. Log er stigende så at maks Like er self også at makse log.
log_lik <- function(rho) {
  temp <- -n/2*log(2*pi) - n/2*log(1-rho^2) - 1/(2*(1-rho^2)) * sum((x_t-rho*x_t_1)^2)
  return(temp)
}

y_s <- sapply(rhos, log_lik)
plot(rhos,y_s) #Hvorfor er den ikke højest i nul så sd er 0.1?


### Nu er det tid til at se den differentierede uden at have ganget med 1-rho^2

rhos2 <- seq(-10.005,9.995, by = 0.01)

log_lik_diff <- function(rho) {
  term1 <- -n/2 * 1/(1-rho^2) * (-2*rho) 
  term2 <- -1/(2*(1-rho^2)) * sum(2*(-x_t_1)*(x_t-rho*x_t_1))
  term3 <- 1/(4*(1-rho^2)^2) * (-4*rho) *sum((x_t-rho*x_t_1)^2)
  return(term1+term2+term3)
}

new_ys <- sapply(rhos2,log_lik_diff)

df <- data.frame(cbind(rhos2,new_ys))

ggplot(df, aes(rhos2, new_ys)) + ylim(-1000,1000) + geom_point()


### Nu er det tid til at se den differentierede EFTER at have ganget med 1-rho^2

rhos2 <- seq(-10.005,9.995, by = 0.01)

log_lik_diff2 <- function(rho) {
  term1 <- (1-rho^2) *n* rho
  term2 <- (1-rho^2) * sum(x_t_1*(x_t-rho*x_t_1))
  term3 <- - rho *sum((x_t-rho*x_t_1)^2)
  return(term1+term2+term3)
}

new_ys <- sapply(rhos2,log_lik_diff2)

df <- data.frame(cbind(rhos2,new_ys))

ggplot(df, aes(rhos2, new_ys)) + ylim(-5000,5000) + geom_point() #Conclusion, den falder hele tiden


### Nu er det tid til at tjekke den differentierede
second_degree_neg <- function(rho) {
  temp <- -2*n*rho^2 + 2*(sum(x_t_1*x_t))*rho+(T-sum(x_t_1^2)-sum(x_t)^2)
  return(temp)
}

new_ys <- sapply(rhos2,second_degree_neg)

df <- data.frame(cbind(rhos2,new_ys))

ggplot(df, aes(rhos2, new_ys)) + geom_point() #Conclusion, den er negativ som forventet.
second_degree_neg(0) #Omkring -8k




