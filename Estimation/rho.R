
est.rho <- function(T_vector,plot = F, arMle=F) {
  #Input: T_vector (as vector or list)
  #Output: rho and m (in a list)
  
  
  #Check input type
  if (typeof(T_vector) == "list") {
    T_vector = unlist(T_vector)
  # } else if (typeof(T_vector) != "vector") {
  } else if (!is.vector(T_vector)) {
    stop("The input of T's has to be either in list or vector")
  }
  
  if(arMle){
    #Estimate rho (Kim's method - standard AR)
    model <- ar.mle(T_vector,aic = F,order.max = 1) #if AIC = F it fits with p = order.max
    
    rho <- model$ar[1]
  } else {
    rho <- optimize(f = est.rhoLoglik, maximum = T,  x = T_vector, lower = -0.99, upper = 0.99)$maximum 
  }
  
  m <- length(T_vector)
  
  #Mangler at add'e den fra AR
  if (plot == T) {
    acf(T_vector,lag.max = max(100,length(T_vector/3)),type = "correlation",plot = T)
  }
  
  return(list(rho = rho,m = m))
    
}

est.rhoLoglik<-function(rho, x){
  N<-length(x)
  
  term1 <- ((N-1)/2)*log(1-rho^2) #Note 0 index convention
  term2 <- (1/(2*(1-rho^2))) * sum((x[2:N]-rho*x[1:(N-1)])^2)
  return(-(term1+term2))
}