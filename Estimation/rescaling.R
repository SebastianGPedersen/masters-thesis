
est.rescale.mu <- function(mu_vector, time_points, t_beginning, bandwidth) {
  
  K2 <- 0.5
  x <- (t_beginning-time_points)/bandwidth
  scaling_integral <- 0.5 * (1-exp(2*x))
  
  new_mu_vector <- sqrt(K2) / sqrt(scaling_integral) * mu_vector
  
  return(new_mu_vector)
}

#Note that sigma2 is sigma^2
est.rescale.sigma <- function(sigma2_vector, time_points, t_beginning, bandwidth) {
  
  K2 <- 0.5
  x <- (t_beginning-time_points)/bandwidth
  scaling_integral <- 0.5 * (1-exp(2*x))
  
  new_sigma2_vector <- K2 / scaling_integral * sigma2_vector
  
  return(new_sigma2_vector)
}