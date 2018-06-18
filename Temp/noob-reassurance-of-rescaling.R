setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("simulation/jumps.R")
source("estimation/estimates.R")
source("estimation/teststat.R")
source("estimation/rho.R")
source("estimation/rescaling.R")
source("kernels/kernels.R")

heston <- sim.heston(sim.setup(Npath = 2))

heston$Y <- t(diff(t(as.matrix(heston$Y))))

hd <- 300/(54*7*24*60*60)
hv <- 5*hd

tind <- seq(20, 23400, 5)


mu<-est.mu.mat.next(heston, hd = hd, t.index = tind)
sig <- est.sigma.mat.next(heston, hv = hv, t.index = tind, lag = 10)

# rescale test

mu_res<-est.rescale.mu(mu$mu, mu$time, 0, hd)
sig_res <- est.rescale.sigma(sig$sig, sig$time, 0, hv)

path <- list(time = heston$time, Y = heston$Y[1,])

mu_vec <- est.mu.next.cpp(path, hd = hd, t.index = tind)
sig_vec <- est.sigma.next.cpp(path, hv = hv, t.index = tind,lag = 10)

mu_res_vec <- est.rescale.mu.vec(mu_vec$mu, mu_vec$time, 0, hd)
sig_res_vec <- est.rescale.sigma.vec(sig_vec$sig, sig_vec$time, 0, hv)

# pre-rescale tests
mu$mu[1,]-mu_vec$mu
sig$sig[1,] - sig_vec$sig

mu_res[1,]-mu_res_vec
sig_res[1,]-sig_res_vec
