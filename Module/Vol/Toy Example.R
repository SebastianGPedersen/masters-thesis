#### source functions neeeded
cd_UNIQUE_NAME_Toy <- getwd()

setwd(paste0(Sys.getenv("masters-thesis"), "\\Vol"))
source("vol_estimators.R")
source("FlexibleFourierForm_Func.R")
source("BV_Analysis_Func.R")
source("PR_Func.R") # PR = Persistence and roughness


setwd(cd_UNIQUE_NAME_Toy)
rm(cd_UNIQUE_NAME_Toy)
##################################

########## DATA
dt <- BV.get_SPY_data(maxDate = "2014-02-01") #Gets first month of data

########### Vol estimation
bucketLengthInMinutes <- 5 # Can't trust results if going lower
bvS_List <- BV.data_deseason_BV_Func(dt = dt, bucketLengthInMinutes = bucketLengthInMinutes)
#(bvS = BV^* in latex)

# Now bvS_list constains results on:
  # Number of terms used in flexible fourier form to deseasonalize data (optimalP)
  # The Buckets used given by their endpoint (e.g. 2014-01-10 09:35) will be 5 minutes are the data starts
  # Number of obervations per bucket
  # A large data table bvSDTfff (fff = flexible fourier form)

print(bvS_List$bvSDTfff) #Quite wide
 
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
  
beta <- PR.est.beta(persistenceDT = persistenceDT, TradingDayLagMin = TradingDayLagMin, TradingDayLagMax = TradingDayLagMax, bucketLengthInMinutes = bucketLengthInMinutes) # OLS is the nicer method, suggested use.
beta

