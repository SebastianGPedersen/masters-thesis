setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")
source("Estimation/estimates.R")
source("estimation/pre-average.R")


#################### PARAMETERS ####################

# Time parameters
mat <- 1
t <- 0:10000  # time_points
dt <- mat/t[length(t)]

# Define sigma
sig <- 1
sig2 <- sig^2
omega <- 0.001
omega2 <- omega^2
ksq <- 0.5 # K2
hd <- 0.1 #bandwidth in mu
hv <- 0.1 #bandwidth in sigma

#Set k_n and phi
theta <- 1/5
k_n <- floor(theta*1/dt^(1/2))
phi_1 <- 1 #int(g'(x)^2)
phi_2 <- 1/12 #int(g^2)
int_g <- 1/4 #int(g)


#################### TESTS ####################

############ WITHOUT NOISE ############

#Sigma
x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
x <- cumsum(x)
dx <- diff(x)

#Pre_avg
pre_x <- est.NewPreAverage(dx,k_n)

#Calculate sigma_hat
C_hat <- 1/(k_n*phi_2)*sum(pre_x^2) #Scaling som forventet
C_hat #Works fine w. 10k obs.


#Mu
N <- 1000
mu_hat <- vector(length = N)

for (i in 1:N) {
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x) 
  dx <- diff(x)
  
  pre_x <- est.NewPreAverage(dx,k_n)
  
  mu_hat[i] <- length(x)/ (length(x)-k_n)*1/(k_n*int_g)*sum(pre_x) #k_n i nævner
}

mean(mu_hat)
var(mu_hat) #Works fine w. 10k obs.


############ WITH NOISE ############

#Sigma
N <- 100 #Takes 3min with N = 100 and t = 100k
sigma_hat <- vector(length = N)

for (i in 1:N) {
  #i <- 1
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  eps <- rnorm(length(t), mean = 0, sd = sqrt(omega^2))
  x <- cumsum(x) 
  y <- x + eps
  dy <- diff(y)
  
  pre_x <- est.NewPreAverage(dy,k_n)
  
  #Calculate sigma_hat
  C_hat <- 1/(k_n*phi_2)*sum(pre_x^2) - phi_1*dt/(2*theta^2*phi_2)*sum(dy^2)
  sigma_hat[i] <- C_hat #Same as Jacod et. al.
}

mean(sigma_hat) #Works fine
var(sigma_hat)

#Mu
N <- 1000 #Takes 3min with N = 100 and t = 100k
mu_hat <- vector(length = N)

for (i in 1:N) {
  
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  eps <- rnorm(length(t), mean = 0, sd = sqrt(omega))
  x <- cumsum(x) 
  y <- x + eps
  dy <- diff(y)
  
  pre_x <- est.NewPreAverage(dy,k_n)
  
  mu_hat[i] <- length(x)/ (length(x)-k_n)*1/(k_n*int_g)*sum(pre_x)
}

mean(mu_hat)
var(mu_hat) #Works fine




