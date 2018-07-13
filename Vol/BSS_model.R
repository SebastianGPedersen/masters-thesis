############################ Base Functions

############## POWER

BSS.g_power <- function(x, alpha, memory_param){
    x^(alpha)*(1+x)^(-(memory_param+alpha))  
}

BSS.power_var <- function(alpha, memory_param){
  beta(2 * alpha + 1, 2 * memory_param - 1)
}

BSS.power_cov <- function(hVec, alpha, memory_param, lower = 0, upper = Inf){
  # sapply to "vectorize" in hVec argument
  sapply(hVec, function(h) {integrate(f = function(x) {BSS.g_power(x+h, alpha, memory_param)*BSS.g_power(x, alpha, memory_param)}, lower = lower, upper = upper, rel.tol = 1e-2)[[1]]})
}

BSS.power_cor <- function(hVec, alpha, memory_param, lower = 0, upper = Inf){
  BSS.power_cov(hVec, alpha, memory_param, lower, upper)/BSS.power_var(alpha, memory_param)
}

## PARALLEL VERSION, requires cluster set-up, e.g.:
# ncl <- max(detectCores()-1, 1)
# cl <- makeCluster(ncl)
# setDefaultCluster(cl = cl) # for ease of calling functions
# clusterExport(cl, "g_power") # cluster environment is generally empty. This moves g_power to all cluster environments
# on.exit(stopCluster(cl)) # stops cluster when cl goes out of scope for any reason, e.g. an error


require(parallel)
BSS.power_cov2 <- function(hVec, alpha, memory_param, lower = 0, upper = Inf){
  # sapply to "vectorize" in hVec argument
  parSapply(cl = NULL, hVec, function(h) {integrate(f = function(x) {BSS.g_power(x+h, alpha, memory_param)*BSS.g_power(x, alpha, memory_param)}, 
                                                    lower = lower, upper = upper, rel.tol = 1e-2)[[1]]})
}

BSS.power_cor2 <- function(hVec, alpha, memory_param, lower = 0, upper = Inf){
  BSS.power_cov2(hVec, alpha, memory_param, lower, upper)/BSS.power_var(alpha, memory_param)
}


############## GAMMA

BSS.g_gamma <- function(x, alpha, memory_param){
    x^(alpha)*exp(-memory_param * x)
}

BSS.gamma_var <- function(alpha, memory_param){
  temp <- 2 * alpha + 1
  gamma(temp)*(2*memory_param)^(-(temp))
}

BSS.gamma_cov <- function(hVec, alpha, memory_param){
  temp <- alpha + 1/2
  gamma(temp+1/2)*(hVec/(2*memory_param))^(temp)*besselK(x = memory_param*hVec, nu = temp)/sqrt(pi)
}

BSS.gamma_cor <- function(hVec, alpha, memory_param){
  BSS.gamma_cov(hVec, alpha, memory_param)/BSS.gamma_var(alpha, memory_param)
}

  
############################ SHARED Functions
# iteration <- 0
optimF <- function(params, hVec, emp, BSS_Cor) {
  alpha <- params[1]
  memory_param <- params[2]
  # iteration <<- iteration + 1
  # print(iteration)
  # print(params)
  if(identical(BSS_Cor, BSS.power_cor)){
    limit <- 0.501
  } else {
    limit <- 0.0000001
  }
  
  if(abs(alpha) >0.499 || memory_param< limit) {
    return(999999)
    
  } else {
    
    res <- tryCatch(
      {sum((BSS_Cor(hVec, alpha, memory_param) - emp)^2)},
      error=function(cond) {return(999999)},
      warning=function(cond) {return(999999) })
    
    return(res)  
    
  }
}


