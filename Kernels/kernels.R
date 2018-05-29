
###############################################
#                                             #
#  A kernel should be specified with as list  #
#  containing kernel function and ksq value   #
#  the ksq value can be NA if not used/known  #
#                                             #
###############################################

# how do we prevent this function from being sourced on its own?
# We could source in source functions, but what happens if we source something that sources?


leftexpkernfunction<-function(x){
  # Input has to be nicely 'sorted' such that lowest nr comes first
  # maybe make kernel a list of function and ksq value?
  neg<-exp(-abs(x[x<=0]))                                                      #only the negatives
  pos<-x[x>0]*0
  return(c(neg,pos))
}

kern.leftexp<-list(kern = leftexpkernfunction, ksq = 0.5)

parzenkernfunction<-function(x, lag){
  #Changed the function. It's extremely annoying that the function changes the order... /Sebastian
  
  x <- x/(lag+1)
  
  x_output <- numeric(length(x))
  
  x_output[x <= 1/2] <- 1-6*x[x <= 1/2]^2+6*x[x <= 1/2]^3
  x_output[x > 1/2] <- 2*(1-x[x > 1/2])^3

  return(x_output)
  
  # # Input has to be nicely 'sorted' such that lowest nr comes first
  # x <- x/(lag+1)
  # xn <- abs(x[x<0])
  # # handles negatives
  # part1n <- 1-6*xn[xn<=0.5]^2+6*xn[xn<=0.5]^3
  # part2n <- 2*(1-xn[xn>0.5 & xn<=1])^3
  # part3n <- xn[xn > 1]*0
  # 
  # # handles positives
  # xp <- abs(x[x>=0])
  # part1p <- 1-6*xp[xp<=0.5]^2+6*xp[xp<=0.5]^3
  # part2p <- 2*(1-xp[xp>0.5 & xp<=1])^3
  # part3p <- xp[xp > 1]*0
  # 
  # return(c(part3n,part2n,part1n,part1p,part2p,part3p))
}

kern.parzen<-list(kern = parzenkernfunction, ksq = NA)