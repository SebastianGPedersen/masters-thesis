setwd(Sys.getenv("masters-thesis"))
source("Kernels/kernels.R")
source("Estimation/estimates.R")
source("estimation/pre-average.R")


#################### PARAMETERS ####################

# Time parameters
mat <- 1
t <- 0:50000  # time_points
dt <- mat/t[length(t)]

# Define sigma
sig <- 1
sig2 <- sig^2
omega <- 0.1
omega2 <- omega^2
ksq <- 0.5 # K2
hd <- 0.1 #bandwidth in mu
hv <- 0.1 #bandwidth in sigma

#Set k_n and phi
k_n <- floor(1/5*1/dt^(1/2))
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
C_hat


#Mu
N <- 100
mu_hat <- vector(length = N)

for (i in 1:N) {
  x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
  x <- cumsum(x) 
  dx <- diff(x)
  
  pre_x <- est.NewPreAverage(dx,k_n)
  
  mu_hat[i] <- length(x)/ (length(x)-k_n)*1/(k_n*int_g)*sum(pre_x) #k_n i nævner
}
#(n-k_n)/n for small sample + epsilon korrection + edge korrection?

mean(mu_hat)
var(mu_hat) # Var = 0.73 w. 500 obs. Var = 0.84 w. 5k obs. Var = 1.38 w. 50k obs. Var = 1.06 w. 100k obs.
#Den stiger m. t. Nok pga. small sample errors

############ WITH NOISE ############

#Sigma
x <- rnorm(n = length(t) , mean = 0, sd = sqrt(sig2*dt))
eps <- rnorm(length(t), mean = 0, sd = sqrt(omega))
x <- cumsum(x)
y <- x + eps

dy <- diff(y)
pre_x <- est.NewPreAverage(dy,k_n) #0.1 sek w. k_n = 100

#Calculate sigma_hat
C_hat <- 1/(k_n*phi_2)*sum(pre_x^2) 
C_hat #1.11 w. 1k


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
mean(mu_hat^2) #0.985 w. 50k after multiplication. k_n = 44
#If n-k_n/n is multiplied, then it decreases w. t 




