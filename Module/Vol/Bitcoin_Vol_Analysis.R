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
##################################

########## DATA
name <- "bitfinex"
# name <- "bitMEX"
# name <- "Kraken"
# maxTol <- 120 # max second without observations
dt <- data.getbitdata(name) #Gets all data
rowsToInclude <- data.Changed(dt$logPrice)
dt <- dt[rowsToInclude]
# dt <- cont.data(DT = dt, maxTol = maxTol) # reduces to longes continous period

########### Vol estimation
bucketLengthInMinutes <- 5 # Can't trust results if going lower than 5
dayLengthInMinutes <- 24*60
bvS_List <- BV.data_deseason_BV_Func(dt = dt, bucketLengthInMinutes = bucketLengthInMinutes, dayLengthInMinutes = dayLengthInMinutes, SPY_bool = F)
#(bvS = BV^* in latex)
min(bvS_List$countDT$N)
# Now bvS_list constains results on:
  # Number of terms used in flexible fourier form to deseasonalize data (optimalP)
  # The Buckets used given by their endpoint (e.g. 2014-01-01 09:35) will be 5 minutes are the data starts
  # Number of obervations per bucket
  # A large data table bvSDTfff (fff = flexible fourier form)

print(bvS_List$bvSDTfff) #Quite wide
plot(bvS_List$bvSDTfff$sCorrect[1:(2*48)], type = "l")
lines(bvS_List$bvSDTfff$nCorrect[1:(2*48)])
plot(bvS_List$bvSDTfff$LogVolsCorrect[1:(2*48)])
plot(bvS_List$bvSDTfff$LogVolnCorrect[1:(2*48)])
# Contains
  # The BV^* estimates in column bvS1 (note the 1 denotes the path, i.e. data viewed as a single path)
  # The volatility (simple calculated found in both theory an BV.data_deseason_BV_Func)
  # Bucket numbers, both day and intraday
  # Log vol
  # corrections: s for simple and n for nice (nice = fff)
    # Described in theory, LogVolsCorrect and LogVolnCorrect are the deseasonalized log vol. 
    # Would always use LogVolnCorrect when needed
  # lagInd used to denote number of lags (here 5 minute periods) if night is not modelled as 0. Use with GREAT caution - it was only created to make a point.



####### Roughness estimation
m <- 6 # Small lag paramter used in estimation. Described in theory. 
variogramDT <- PR.variogram_prep_DT(bvS_List$bvSDTfff) # Extract (and rename) columns

alpha <- PR.est.alpha(variogramDT = variogramDT, m = m, bucketLengthInMinutes = bucketLengthInMinutes, OLS = T) # OLS is the nicer method, suggested use.
alpha

####### Persistence estimation
persistenceDT <- PR.persistence_prep_DT(bvS_List$bvSDTfff) # Extract (and rename) columns
TradingDayLagMin <- floor(nrow(persistenceDT)^(1/4)) # Described in theory. Not the suggested value, just used for demonstration
TradingDayLagMax <- floor(nrow(persistenceDT)^(1/3)) # Described in theory. Not the suggested value, just used for demonstration
  
beta <- PR.est.beta(persistenceDT = persistenceDT, TradingDayLagMin = TradingDayLagMin, TradingDayLagMax = TradingDayLagMax, bucketLengthInMinutes = bucketLengthInMinutes, SPY_Bool = F) # OLS is the nicer method, suggested use.
beta

