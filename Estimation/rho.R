#The first function is MLE and the one we should use
#The second is Kim's OLS estimate

est.rho.MLE <-function(T_vector, error = 0.0001){ #With newest observation at index 1
  
  #T_vector <- c(1,2,3,2,4,5,1,2)
  X_t <- T_vector[-length(T_vector)] 
  X_t_1 <- T_vector[-1]
  
  N <- length(X_t)
  
  #Polynomium
  third_order_pol <- function(p) {
    value <- (1-p^2)*N*p + (1-p^2)*sum(X_t_1 * (X_t - p*X_t_1)) - p * sum((X_t - p*X_t_1)^2)
    return(value)
  }
  
  #Bisection (with p = -1 the function is positive, with p = 1 it is negative)
  
  lower <- -1
  upper <- 1
  
  while (upper-lower > 2*error) {
    mid <- (upper+lower)/2
    temp_val <- third_order_pol(mid)
    if (temp_val > 0) {lower <- mid} else {upper <- mid}
  }

  rho_est <- (upper+lower)/2
  
  return(rho_est)
}

est.rho.kim <- function(T_vector,plot = F, arMle=F) {
  #Input: T_vector (as vector or list)
  #Output: rho and m (in a list)
  
  #Check input type
  if (typeof(T_vector) == "list") {
    T_vector = unlist(T_vector)
    # } else if (typeof(T_vector) != "vector") {
  } else if (!is.vector(T_vector)) {
    stop("The input of T's has to be either in list or vector")
  }
  
  #Estimate rho (Kim's method - standard AR)
  model <- ar.mle(T_vector,aic = F,order.max = 1) #if AIC = F it fits with p = order.max
  
  rho <- model$ar[1]
  sigma <- model$var.pred
  m <- length(T_vector)
  
  #Mangler at add'e den fra AR
  if (plot == T) {
    acf(T_vector,lag.max = max(100,length(T_vector/3)),type = "correlation",plot = T)
  }
  
  return(rho)
  
}

est.rho.third_degree <- function(T_vector) {
  
  ### Calculate the three roots
  N <- length(T_vector)
  x_t_1 <- T_vector[1:(N-1)]
  x_t <- T_vector[2:N]
  
  a <- -N
  b <- sum(x_t*x_t_1)
  c <- N - sum(x_t_1^2) - sum(x_t^2)
  d <- b
  
  
  Delta_0 <- as.complex(b^2-3*a*c)
  Delta_1 <- 2*b^3 - 9*a*b*c + 27*a^2*d
  C_1 <- ((Delta_1 + sqrt(Delta_1^2-4*Delta_0^3))/2)^(1/3)
  
  xi <- (-1/2+1/2*sqrt(3)*complex(imaginary = 1))
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
    for (i in all_x) {
      if (abs(Im(i)) < 10^(-10)) {return(Re(i))}
    }
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

