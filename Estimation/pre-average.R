est.PreAverage<-function(logRet_Data, K){
  if((!(K%%2==0))&(K>1)){
    stop("Needs K even")
  }
  
  N <- length(logRet_Data) 
  lOut <- N-K+1
  r <- rep(NA, lOut)
  
  for(i in 1:(lOut)){
    
    r[i] <- 1/K * (sum(logRet_Data[i+(K/2):(K-1)])-sum(logRet_Data[i+(0):(K/2-1)]))
  }  
  
  #Zeros at first
  r_w_zeros <- c(rep(0,K-1),r)
  
  return(r_w_zeros/K)
}

# Consider adding k-1 zeros in the beginning


#delta_y <- dy[,1]
#k_n <- k_n_list[1]
#kernel_func = parzen_kernel_func
#The function below has to be quick, so the loop is in counter-intuitive order (k_n is outer and n is inner loop)
est.NewPreAverage <- function(delta_y, k_n, kernel_func = NA) {

  #If kernel function is missing, let it be standard
  if (is.na(kernel_func)) {
    kernel_func <- function(x) {min(x,1-x)}
  }
  
  #Calculate kernel_values outside loop
  kernel_values <- sapply((1:(k_n-1))/k_n, kernel_func)
  
  #Initialize
  y_bar <- rep(0,length(delta_y)-k_n+2)
  
  #Loop over k_n
  for (i in 1:(k_n-1)){
    #i <- 1
    y_bar <- y_bar + kernel_values[i]*delta_y[(k_n-i):(length(delta_y)+1-i)] #længde n-k_n+2
  }

  return(y_bar)
}

