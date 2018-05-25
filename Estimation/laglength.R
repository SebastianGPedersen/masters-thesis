# LAG LENGTH #
# find Newey-West (1994) -- automatic lag selection with parzen kernel

laglength = function(dx, nmu){
  # dx are NOT preavr
  # only dx from backwards in time should be considered
  # Indices may be wrong due to zero index (or is it? Matlab is index 1...)
  c <- 2.6614
  q <- 2
  
  n <- round(4*(nmu/100)^(4/25))
  root <- 1/(2*q+1)
  
  # pre-estimate autocovariance of noise increment
  E_deij <- numeric(n+1)
  for(i in 1:(n+1)){
    #i <- 2
    E_deij[i] <- mean(dx[(i+1):length(dx)]*dx[1:(length(dx)-i)])
  }
  
  v0 <-  -sum((1:(n+1))*E_deij)
  
  ac<- numeric(n)
  
  for(i in 1:n){
    ac[i] <- -sum(   (1:(n+1-i))*E_deij[(i+1):(n+1)]   )
  }
  
  # determine lag length
  
  s0 <- v0 + 2 * sum(ac)
  sq <- 2*sum(   (1:n)^(q)*ac   ) # sq = 2*sum((1:n).^q.*ac); - regnehieraki
  
  gamma <- c*(   ((sq/s0)^q)^root   )
  return(round(gamma*nmu^root))
}