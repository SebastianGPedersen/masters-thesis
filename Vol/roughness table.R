
#### source functions neeeded
cd_UNIQUE_NAME_Toy <- getwd()

setwd(paste0(Sys.getenv("masters-thesis"), "\\Vol"))
source("vol_estimators.R")
source("FlexibleFourierForm_Func.R")
source("BV_Analysis_Func.R")
source("PR_Func.R") # PR = Persistence and roughness
source("adf_fit_func.R") 

setwd(paste0(Sys.getenv("masters-thesis"), "\\SPY"))
source("datahandling.R")

setwd(cd_UNIQUE_NAME_Toy)
rm(cd_UNIQUE_NAME_Toy)


############ DATA
name <- "bitfinex"
# name <- "bitMEX"
# name <- "Kraken"
# maxTol <- 120 # max second without observations
dt <- data.getbitdata(name) #Gets all data
rowsToInclude <- data.Changed(dt$logPrice)
dt <- dt[rowsToInclude]
###########





dayLengthInMinutes <- 24*60
bucketLengthInMinutes <- c(5,10,15,30)
m <- 3:200
NLS_alphaMatrix <- matrix(NA, nrow = length(m), ncol = length(bucketLengthInMinutes)) 
OLS_alphaMatrix <- matrix(NA, nrow = length(m), ncol = length(bucketLengthInMinutes)) 
PVec <- rep(NA, length(bucketLengthInMinutes))

for(i in seq_along(bucketLengthInMinutes)){
  print(i)
  
  #Compute IV given bucket lengths
  temp <- BV.data_deseason_BV_Func(dt = dt, bucketLengthInMinutes = bucketLengthInMinutes[i], dayLengthInMinutes = dayLengthInMinutes, SPY_bool = F)
  bvSDTfff <- temp$bvSDTfff
  PVec[i] <- temp$optimalP

  rm(temp)

  #Quick extract/rename columns
  variogramDT <- PR.variogram_prep_DT(bvSDTfff)

  
  for(j in seq_along(m)){
    print(c(i, j))
    #Compute alpha
    NLS_alphaMatrix[j, i] <- PR.est.alpha(variogramDT = variogramDT, m = m[j], bucketLengthInMinutes = bucketLengthInMinutes[i], OLS = F, SPY_Bool = F)
    OLS_alphaMatrix[j, i] <- PR.est.alpha(variogramDT = variogramDT, m = m[j], bucketLengthInMinutes = bucketLengthInMinutes[i], OLS = T, SPY_Bool = T)
  }
}

# NLS_alphaMatrix <- readRDS(paste0("NLS_AlphaTable__Buckets_",
#                                                        min(bucketLengthInMinutes), "-", max(bucketLengthInMinutes),
#                                                        "__Lags_",
#                                                        min(m), "-", max(m),
#                                                        ".rds"))
# OLS_alphaMatrix <- readRDS(paste0("OLS_AlphaTable__Buckets_",
#                                   min(bucketLengthInMinutes), "-", max(bucketLengthInMinutes),
#                                   "__Lags_",
#                                   min(m), "-", max(m),
#                                   ".rds"))

# library(plotly)
# p <- plot_ly(z = ~OLS_alphaMatrix[,4:7], x = 1:4 ) %>% add_surface()
# p

library(ggplot2)
library(reshape2)
# prows <- 1:10
prows <- 1:length(m)
# columns <- 1:7
columns <- 1:(length(bucketLengthInMinutes))
cnames <- paste0(bucketLengthInMinutes[columns], " minutes")
nlsdf <- data.frame(NLS_alphaMatrix[prows,columns])
names(nlsdf) <- cnames
nlsdf$Lag <- m[prows]
nlsdf$Type <- "NLS"

nlspdf <- melt(nlsdf, id = c("Lag", "Type"), variable.name = "BucketWidth", value.name = "RoughnessIndex")


olsdf <- data.frame(OLS_alphaMatrix[prows,columns])
names(olsdf) <- cnames
olsdf$Lag <- m[prows]
olsdf$Type <- "OLS"

olspdf <- melt(olsdf, id = c("Lag", "Type"), variable.name = "BucketWidth", value.name = "RoughnessIndex")

pdf <- rbind(nlspdf,olspdf)

# cols1f <- colorRampPalette(c("red", "orange"))
# cols1 <- cols1f(4)
# cols2f <- colorRampPalette(c("blue", "purple"))
# cols2 <- cols2f(4)
# cols <- c(cols1,cols2)
# cols <- cols1f(8)
ggplot(data=pdf,
  aes(x=Lag, y=RoughnessIndex, colour=BucketWidth, shape = Type,
      group=interaction(BucketWidth, Type))) +
  geom_line() +
  geom_point(size = 3) + 
  # scale_color_manual(values = cols) +
  labs(y = "Rougness Index",  color = "Bucket", x = "m") + 
  scale_x_continuous(breaks = m[floor(seq(min(prows),max(prows), length.out = 30))])



