T_vector <- rnorm(n =1000, mean = 0,sd = 0.9)
x_t <- X_es[2:1000]
x_t_1 <- X_es[1:999]
n <- length(x_t)

est.rho.third_degree(T_vector) #Denne siger én løsning

est.rho.third_degree <- function(T_vector) {
  
  ### Calculate the three roots
  N <- length(T_vector)
  x_t_1 <- T_vector[1:(N-1)]
  x_t <- T_vector[2:N]
  
  a <- -N
  b <- sum(x_t*x_t_1)
  c <- N - sum(x_t_1^2) - sum(x_t^2)
  d <- b
  
  
  Delta <- as.complex(18*a*b*c*d - 4*b^3*d + b^2*c^2 - 4*a*c^3 - 27*a^2*d^2)
  Delta_0 <- as.complex(b^2-3*a*c)
  Delta_1 <- 2*b^3 - 9*a*b*c + 27*a^2*d
  
  C_1 <- ((Delta_1 + sqrt(Delta_1^2-4*Delta_0^3))/2)^(1/3)

  x_1 <- -1/(3*a) * (b+xi^0*C_1+Delta_0 /(xi^0*C_1))
  x_2 <- -1/(3*a) * (b+xi^1*C_1+Delta_0 /(xi^1*C_1))
  x_3 <- -1/(3*a) * (b+xi^2*C_1+Delta_0 /(xi^2*C_1))
  
  #order the solutions
  all_x <- c(x_1,x_2,x_3)
  
  #Count number of real solutions
  n_real <- 0
  for (i in all_x) {
    if (abs(Im(i)) < 10^(-10)) {n_real <- n_real+1} #numerical issues
  }
  
  if (n_real == 0 || n_real == 2) {stop("The third degree solution does not work")}
  
  ## If a single solution
  if (n_real == 1) {
    if (abs(Im(i)) < 10^(-10)) {return(Re(x_1))} else if (abs(Im(i)) < 10^(-10)) {return(Re(x_2))} else {return(Re(x_3))}
  }
  
  ### If 3 solutions
  #order solutions
  all_real_x <- c(Re(x_1), Re(x_2), Re(x_3))
  all_real_x <- sort(all_real_x)
  
  
  #Find out if first or last observation maximizes 
  
  log_lik <- function(rho) {
    term1 <- -N/2 * log(2*pi)
    term2 <- -N/2 * log(1-rho^2)
    term3 <- -1/(2*1-rho^2) * sum((x_t-rho*x_t_1)^2)
    return(term1+term2+term3)
  }
  
  log_lik1 <- log_lik(all_real_x[1])
  log_lik2 <-log_lik(all_real_x[3])
  
  if (log_lik1 > log_lik2) {
    return(all_real_x[1])
  } else{
    return(all_real_x[3])
  }
}
