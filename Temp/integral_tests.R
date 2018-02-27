
#Integralet på s. 16 er givet ved 

c_1 <- 3
alpha <- 0.75

mat <- 6.5/(24*7*52)
begin <- 0.475*mat
tau <- 0.5*mat
begin2 <- 0
tau2 <- 0.025*mat

######## SYMBOLIC INTEGRAL (with substitution error)
integral_func <- function(c_1,alpha,tau,begin) {
return <- -c_1 *1/(-alpha+1) * tau^(-alpha + 1) +c_1 *1/(-alpha+1) * (begin)^(-alpha + 1) 
return(return)
}

integral_func(c_1,alpha,tau,begin) #Dette giver 2 %


######## APPROXIMATION

steps <- 10000
dt <- (tau2-begin2)/steps

mu = 0
t = begin2
for (i in 1:steps) {
  mu <- mu - c_1/(tau-t)^alpha * dt
  t <- t + dt
}
mu

# Sustitution error og samme antal steps som Kim

steps <- 23400*(0.5-0.475)
dt <- (tau2-begin2)/steps

mu = 0
t = begin2

for (i in 1:steps) {
  mu <- mu - c_1/(tau-t)^alpha * dt
  t <- t + dt
}
mu

######## TEST
-c_1/(1-alpha)*0.025^(1-alpha) #Dette giver -4.77 (bruger relative tal)
-c_1/(1-alpha)*(0.025*mat)^(1-alpha) #Dette giver -0.78 (bruger korrekte tal)

-c_1/(1-alpha)*(0.5*mat)^(1-alpha)+c_1/(1-alpha)*(0.475*mat)^(1-alpha) #Dette giver -2% (laver substitutionsfejl)



######## TEST AF VOLATILITY INCREASE IN WINDOW

tau <- 0.5*mat
begin <- 0.475*mat
end <- 0.552*mat
begin2 <- 0
end2 <- 0.05

beta <- 0.4
c_2 <- 0.15

### SYMBOLIC INTEGRAL

## correct
upper <- c_2/(1-beta) * end^(1-beta)
lower <- c_2/(1-beta) * begin^(1-beta)
my_int <- upper - lower

(avg_increase <- my_int/(end-begin))

## substitution error
upper <- c_2/(1-beta) * end2^(1-beta)
lower <- c_2/(1-beta) * begin2^(1-beta)
my_int <- upper - lower

(avg_increase <- my_int/(end2-begin2))


## relative numbers error
upper <- c_2/(1-beta) * 0.525^(1-beta)
lower <- c_2/(1-beta) * 0.475^(1-beta)
my_int <- upper - lower

(avg_increase <- my_int/(end2-begin2))
