
leftexpkern<-function(x){
  # Input has to be nicely 'sorted' such that lowest nr comes first
  # maybe make kernel a list of function and ksq value?
  neg<-exp(-abs(x[x<=0]))                                                      #only the negatives
  pos<-x[x>0]*0
  return(c(neg,pos))
}

parzenkern<-function(x){
  # Input has to be nicely 'sorted' such that lowest nr comes first
  xn <- abs(x[x<0])
  # handles negatives
  part1n <- 1-6*xn[xn<=0.5]^2+6*xn[xn<=0.5]^3
  part2n <- 2*(1-xn[xn>0.5 & xn<=1])^3
  part3n <- xn[xn > 1]*0
  
  # handles positives
  xp <- abs(x[x>=0])
  part1p <- 1-6*xp[xp<=0.5]^2+6*xp[xp<=0.5]^3
  part2p <- 2*(1-xp[xp>0.5 & xp<=1])^3
  part3p <- xp[xp > 1]*0
  
  return(c(part3n,part2n,part1n,part1p,part2p,part3p))
}