
est.rho <- function(T_vector,plot = F) {
  #Input: T_vector (as vector or list)
  #Output: rho and m (in a list)
  
  
  #Check input type
  if (typeof(T_vector) == "list") {
    T_vector = unlist(T_vector)
  } else if (typeof(T_vector) != "vector") {
    stop("The input of T's has to be either in list or vector")
  }
  
  #Estimate rho (Kim's method - standard AR)
  model <- ar.mle(T_vector,aic = F,order.max = 1) #if AIC = F it fits with p = order.max
  
  rho <- model$ar[1]
  m <- length(T_vector)
  
  #Mangler at add'e den fra AR
  if (plot == T) {
    acf(T_vector,lag.max = max(100,length(T_vector/3)),type = "correlation",plot = T)
  }
  
  return(list(rho = rho,m = m))
    
}