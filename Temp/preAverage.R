est.preAverage<-function(data, t.index, freq, overlap = F, width = 0){
  if(missing(t.index)){
    
  }
  else{
    
  }
}

data <- list(time = 1:100, Y = rnorm(100, 0, 1))

freq = 30

overlap = F

times <- floor(length(data$time)/freq)
indices <- seq(from = length(data$time)-(times*freq), to = length(data$time) , by=freq )

t.index <- c(10, 30, 50, 90)

start <- fuuuck
