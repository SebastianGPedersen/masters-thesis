setwd(Sys.getenv("masters-thesis"))
source("simulation/heston.R")
source("simulation/bursts.R")
source("estimation/estimates.R")
source("estimation/rho.R")
source("estimation/teststat.R")
source("estimation/pre-average.R")
source("kernels/kernels.R")
source("simulation/jumps.R")

dbVjump <- function(Nstep, burstsetting, index, seed){
  burstsim <- function(settings, burstsetting){
    alpha <- burstsetting$alpha
    beta <- burstsetting$beta
    burst_time <- burstsetting$burst_time
    interval_length <- burstsetting$interval_length
    c_1 <- burstsetting$c_1
    c_2 <- burstsetting$c_2
    
    sims<-sim.heston(setting)
    
    sims.vb<-sim.addvb(sims,    burst_time = burst_time, interval_length = interval_length,
                       c_2 = c_2, beta  = beta)
    
    sims.db<-sim.adddb(sims.vb, burst_time = burst_time, interval_length = interval_length,
                       c_1 = c_1,  alpha = alpha)
    sims.jump <- sim.addjump(sims, burst_time = burst_time, interval_length = interval_length, c_1 = c_1, alpha = alpha)
    
    
    #Get a single path
    path = 1
    
    Heston_path = sim.path(path,sims)$Y
    vb_path = sim.path(path, sims.vb)$Y
    vbdb_path = sim.path(path, sims.db)$Y
    
    
    return(list(raw = sim.path(path,sims), vb = sim.path(path, sims.vb), vbdb = sim.path(path, sims.db), jump = sim.path(path,sims.jump)))
  }
  
  #Sim
  if(!missing(seed)){
    set.seed(seed)
  }
  setting <- sim.setup(Npath = 2, Nstep = Nstep, omega = 0.0000225)
  sims<-burstsim(setting, burstsetting)
  
  
  #prev
  theta <- 1
  k <- theta*( floor(sqrt(length(sims$raw$time)))-floor(sqrt(length(sims$raw$time)))%%2 ) #force to even integer
  prev.hest <- c(rep(0,k-2+2),est.NewPreAverage(diff(sims$raw$Y), k))
  prev.db <- c(rep(0,k-2+2),est.NewPreAverage(diff(sims$vbdb$Y), k))
  prev.jump <- c(rep(0,k-2+2),est.NewPreAverage(diff(sims$jump$Y), k))
  
  data.hest <- list(time = sims$raw$time*1*52*7*24*60*60*1000, Y = prev.hest, raw = sims$raw$Y)
  data.db <- list(time = sims$vbdb$time*1*52*7*24*60*60*1000, Y = prev.db, raw = sims$vbdb$Y)
  data.jump <- list(time = sims$jump$time*1*52*7*24*60*60*1000, Y = prev.jump, raw = sims$jump$Y)
  
  #Estimation
  if(missing(index)){
    tind <- seq(from = 2000, to = Nstep-1, by = 60)
  }
  else{
    tind <- index
  }
  dt <- diff(data.hest$time)[1]
  hd <- 300*dt
  hv <- 300*dt
  #hd<- 3000*sqrt(dt)
  #hv<- 3000*sqrt(dt)
  
  # Calc T
  test.hest <- test.db(data = data.hest, hd = hd, hv = hv, t.index = tind, noisefun = est.noise.iid.next, theta = theta, kn = k)
  test.burst <- test.db(data = data.db, hd = hd, hv = hv, t.index = tind, noisefun = est.noise.iid.next, theta = theta, kn = k)
  test.jump <- test.db(data = data.jump, hd = hd, hv = hv, t.index = tind, noisefun = est.noise.iid.next, theta = theta, kn = k)
  
  # Calc Tstar
  Tstar.hest<-tstar(test.hest)$tstar
  Tstar.burst<-tstar(test.burst)$tstar
  Tstar.jump<-tstar(test.jump)$tstar
  
  return(c(Tstar.hest,Tstar.burst, Tstar.jump))
  
  # fit rho
  #rho.hest <- est.rho(test.hest$test)
  #rho.burst <- est.rho(test.burst$test)
  #rho.jump <- est.rho(test.jump$test)
  
  #
  #conf <- 0.95
  #z.hest<-est.z_quantile(rho.hest$m, rho.hest$rho, conf)$qZm
  #z.burst<-est.z_quantile(rho.burst$m, rho.burst$rho, conf)$qZm
  #z.jump<-est.z_quantile(rho.jump$m, rho.jump$rho, conf)$qZm
  
  #res<-c(Tstar.hest>=z.hest, Tstar.burst>=z.burst, Tstar.jump>=z.jump)
  #names(res) <- c("hest", "burst", "jump")
}

# Test N's
N <- c(5.6, 10.4, 20.4, 30, 40.4, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190,
       200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330.4, 340.6, 350.6)*1000

# setup matrix
res <- matrix(NA, length(N), 3)
colnames(res) <- c("hest", "burst", "jump")

res2 <- matrix(NA, length(N), 3)
colnames(res2) <- c("hest", "burst", "jump")

# settings
burstset<-sim.burstsetting(alpha = 0.65, beta = 0.1 ,c_1 = 0.1, c_2 = 0.02, interval_length = 0.1)
burstset2<-sim.burstsetting(alpha = 0.7, beta = 0.1 ,c_1 = 0.1, c_2 = 0.02, interval_length = 0.1)

for(i in 1:length(N)){
  try(res[i,] <- dbVjump(N[i], burstset, seed = 12444, index = N[1]/2+floor(sqrt(N[1])) ))
  
  try(res2[i,] <- dbVjump(N[i], burstset2, seed = 12444,index = N[1]/2+floor(sqrt(N[1])) ))
  print(i)
}

res<-cbind(N,res)
res2<-cbind(N,res2)

clean.res<-res[complete.cases(res),]
clean.res2<-res2[complete.cases(res2),]

require(ggplot2)
data1 <- data.frame(clean.res)
data2 <- data.frame(clean.res2)

ggplot() +
  geom_line(data=data1, aes(x=N, y=hest, col = "hest - small"), size = 1) +
  geom_line(data=data1, aes(x=N, y=burst, col = "burst - small"), size = 1) +
  geom_line(data=data1, aes(x=N, y=jump, col = "jump - small"), size = 1) +
  geom_line(data=data2, aes(x=N, y=hest, col = "hest - large"), size = 1) +
  geom_line(data=data2, aes(x=N, y=burst, col = "burst - large"), size = 1) +
  geom_line(data=data2, aes(x=N, y=jump, col = "jump - large"), size = 1) +
  xlab("Nsteps") + ylab("T*")  + 
  ggtitle("T* as function of Nsteps")

#test if Nsteps is compatible with sim burst
#dbVjump(350600, burstset, seed = 12444)

if(0 > 1){
  saveRDS(res, file=paste0("temp/","hd=dt^0.25","alpha=",burstset$alpha,"beta=",burstset$alpha,"c1=",burstset$c_1,"c2=",burstset$c_2,".Rda"))
  saveRDS(res2, file=paste0("hd=dt^0.25","alpha=",burstset2$alpha,"beta=",burstset2$alpha,"c1=",burstset2$c_1,"c2=",burstset2$c_2,".Rda"))
}

  paste0("temp/","hd=","alpha=",burstset$alpha,"beta=",burstset$alpha,"c1=",burstset$c_1,"c2=",burstset$c_2,".Rda")
  paste0("hd=","alpha=",burstset2$alpha,"beta=",burstset2$alpha,"c1=",burstset2$c_1,"c2=",burstset2$c_2,".Rda")
  
# Hard notify when done
#shell.exec("https://www.youtube.com/embed/quxTnEEETbo?autoplay=1")

