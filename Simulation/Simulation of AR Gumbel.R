# SIMULATION OF QUANTILE (AR GUMBEL)
setwd(Sys.getenv("masters-thesis-data"))

rhoset <- 0:9/10
mset <- seq(5, 100, 10)
n <- 10000
Npath = 1*10^6
q<- 0.95

# FUNCTION START
star <- function(x) abs(max(x))
# 
out <- matrix(NA, nrow = length(rhoset), ncol = length(mset))
inner <- matrix(NA, nrow = Npath, ncol = length(mset))
# for rho starts #
for(r in 1:length(rhoset)){
  rho <- rhoset[r]
  print(rho)
  # for each path #
  for(i in 1:Npath){
    X<-arima.sim(model=list(ar=rho), n = n)
    end <- length(X)
    
    for(j in 1:length(mset)){
       inner[i,j]<- star(X[(end-mset[j]):end])
    }
  }
  for(j in 1:length(mset)){
    out[r,j]<-as.numeric(quantile(inner[,j], probs = q))
  }
  save(out, file = "smallm.RData")
}


if(F){
  ja <- load("smallm.RData")
  
  assign("interpolList", readRDS("interpolList.rds"), envir = .GlobalEnv)
  
}

interpol <- readRDS("interpolList.rds")
