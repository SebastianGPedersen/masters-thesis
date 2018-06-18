#time_points <- Heston$time[desired_indices]
#t_beginning <- Heston$time[1]
#h_sigma <- h_mu*ratio

est.rescale.mu <- function(mu_matrix, time_points, t_beginning, h_mu) {
  
  K2 <- 0.5
  x <- (t_beginning-time_points)/h_mu
  scaling_integral <- 0.5 * (1-exp(2*x))
  
  new_mu_matrix <- mu_matrix
  
  #This is the fastest wa
  for (col in 1:ncol(mu_matrix)) {
    new_mu_matrix[,col] <- mu_matrix[,col] / sqrt(scaling_integral[col])
  }
  
  new_mu_matrix <- sqrt(K2) * new_mu_matrix
  return(new_mu_matrix)
}

#Note that sigma2 is sigma^2
est.rescale.sigma <- function(sigma2_matrix, time_points, t_beginning, h_sigma) {
  
  K2 <- 0.5
  x_2 <- (t_beginning-time_points)/h_sigma
  scaling_integral_2 <- 0.5 * (1-exp(2*x_2))
  
  new_sigma2_matrix <- sigma2_matrix
  
  for (col in 1:ncol(sigma2_matrix)) {
    new_sigma2_matrix[,col] <- K2 / scaling_integral_2[col] * sigma2_matrix[,col]
  }
  
  return(new_sigma2_matrix)
}

# VECTOR BASED FOR REAL DATA
est.rescale.mu.vec <- function(mu_vec, time_points, t_beginning, h_mu) {
  
  K2 <- 0.5
  x <- (t_beginning-time_points)/h_mu
  scaling_integral <- 0.5 * (1-exp(2*x))
  
  #This is the fastest wa
  new_mu_vec <- mu_vec/ sqrt(scaling_integral)
  
  new_mu_vec <- sqrt(K2) * new_mu_vec
  return(new_mu_vec)
}

#Note that sigma2 is sigma^2
est.rescale.sigma.vec <- function(sigma2_vec, time_points, t_beginning, h_sigma) {
  
  K2 <- 0.5
  x_2 <- (t_beginning-time_points)/h_sigma
  scaling_integral_2 <- 0.5 * (1-exp(2*x_2))
  
  new_sigma2_vec <- K2 / scaling_integral_2 * sigma2_vec
  
  return(new_sigma2_vec)
}