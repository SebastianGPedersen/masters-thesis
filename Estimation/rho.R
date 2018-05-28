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
  
  T_vector <- X
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
  
  return(list(rho = rho,sigma = sigma))
  
}
