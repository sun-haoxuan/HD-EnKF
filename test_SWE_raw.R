rm(list = ls())
setwd("E:/Project/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)

source('Code/Functions/cal_GTE_bandwidth.R')
source('Code/Functions/cal_HPHR.R')
source('Code/Functions/cal_RMSE.R')
source('Code/Functions/ev_adjust.R')
source('Code/Functions/SWE_State.R')
source('Code/Functions/SWE_Forecast.R')
source('Code/Functions/SWE_Analyse.R')
source('Code/Functions/RK4.R')
source('Code/Functions/taper_functions.R')
source('Code/Functions/dist_functions.R')
source('Code/cal_taper_cov_SWE.R')
# Rcpp::sourceCpp('Code/Functions/cal_bandwidth_m0.cpp')

set.seed(1234)
opt = list(
  ## Generate true state
  nx = 50, # grid size of x
  ny = 31, # grid size of ys
  L = 5e+5,
  D = 3e+5,
  g = 9.8, # gravitational constant
  f = 1e-4,
  k = NA,
  dt = 30, # hardwired timestep
  h0 = 50,
  h1 = 5.5,
  h2 = 3.325,
  ## Generate observations
  # id = c(5,9,16,21,25,30,31,33,44,49),
  id = sort(sample(50, 10)),
  sigma.obs = c(sqrt(0.5), sqrt(0.5), 1), # observation error
  rho = 0.5, ## observation
  vf = 1, # for Wrong R
  ## Assimilation Step
  S.o.seq = 360, # assimilate every 3h
  S.o.start = 61, ## Start of the observation
  S.o.end = 5820, # assimilate 48h, then forecast 24h
  S = 8700,
  ## Assimilation
  n = 100,
  k.true = 1e+4,
  k.set = 5e+4,
  sigma.i = 1,
  method = 'standard'
)
state = SWE_State(opt)
analyse = SWE_Analyse(state, opt)

t1 = Sys.time()

opt$k = opt$k.set

## Initialize ####
x.t = state$true
y.o = state$observation

## Generate true state
nx = opt$nx # grid size of x
ny = opt$ny # grid size of y
L = opt$L
D = opt$D
g = opt$g # gravitational constant
f = opt$f
k = opt$k # given by function
dt = opt$dt # hardwired timestep
h0 = opt$h0
h1 = opt$h1
h2 = opt$h2

## Generate observations
id = opt$id # rows with observation
sigma.o = opt$sigma.o # observation error
rho = opt$rho ## observation
vf = opt$vf # for Wrong R

## Assimilation Step
S = opt$S
S.o.start = opt$S.o.start ## Start of the observation
S.o.end = opt$S.o.end # assimilate 48h, then forecast 24h
S.o.seq = opt$S.o.seq # assimilate every 3h

## Assimilation
n = opt$n
k.set = opt$k.set
sigma.i = opt$sigma.i
method = opt$method

## Other parameters
ngrid <- nx * ny
dx <- L / nx
dy <- D / (ny - 1)
nobs <- (S.o.end - S.o.start) %/% S.o.seq + 1 # forecast 1h, then assimilate
ndim <- ngrid * 3
nid <- length(id)
xobs <- ny * nid * 3 # observe h only
coor = data.frame(
  x = rep(1:nx, times = ny),
  y = rep(1:ny, each = nx)
)
ncoor = nrow(coor)
digits = 1
mat.dist = diag(0, ncoor)
for(i in 1:(ncoor - 1)){
  for(j in (i + 1):ncoor){
    d_ij = round(sqrt(sum((coor[i, ] - coor[j, ]) ^ 2)), digits)
    mat.dist[i, j] = mat.dist[j, i] = d_ij
  }
}

## Observation Matrix ####
H <- matrix(0, nrow = xobs, ncol = ndim)
for (i in 1:nid) {
  jy <- (id[i] - 1) * ny + 1:ny
  jx <- ((i - 1) * ny + 1):(i * ny)
  for (ii in 1:ny) {
    H[jx[ii], jy[ii]] <- 1
    H[jx[ii] + ny * nid, jy[ii] + ngrid] <- 1
    H[jx[ii] + ny * nid * 2, jy[ii] + ngrid * 2] <- 1
  }
}
rm(jx, jy)

## Observation error covariance matrix ####
aa1 <- matrix(NA, nrow = 2, ncol = ny * nid)
for (i in 1:nid){
  iy <- ((i - 1) * ny + 1):(i * ny)
  for (ii in 1:ny) {
    aa1[1, iy[ii]] <- ii
    aa1[2, iy[ii]] <- id[i]
  }
}
aa2 <- diag(1, ny * nid)
for (i in 1:(ny * nid - 1)){
  for (j in (i + 1):(ny * nid)) {
    dd <- sqrt((aa1[1, i] - aa1[1, j])^2 + (aa1[2, i] - aa1[2, j])^2)
    aa2[i, j] <- rho^dd
    aa2[j, i] <- aa2[i, j]
  }
}

R <- matrix(0, nrow = xobs, ncol = xobs)
R[1:(ny * nid), 1:(ny * nid)] <- aa2 * vf * sigma.o[1]^2
R[(ny * nid + 1):(ny * nid * 2), (ny * nid + 1):(ny * nid * 2)] <- aa2 * vf * sigma.o[2]^2
R[(ny * nid * 2 + 1):(ny * nid * 3), (ny * nid * 2 + 1):(ny * nid * 3)] <- aa2 * vf * sigma.o[3]^2

R.svd <- svd(R)
R.inv <- solve(R)
R.inv.root = R.svd$v %*% diag(R.svd$d ^ (- 0.5)) %*% t(R.svd$u)
logdetR <- sum(log(R.svd$d))

rm(aa1, aa2, R.svd)

## Observation State ####
y.f.en <- array(NA, c(xobs, n, nobs))
t.o <- rep(0, times = S)
for (i in seq(S.o.start, S.o.end, by = S.o.seq)) {
  t.o[i] <- 1
}
so <- H %*% x.t[, t.o == 1]

R.svd.u <- svd(R[1:(ny * nid), 1:(ny * nid)] / vf)
R.svd.v <- svd(R[(ny * nid + 1):(ny * nid * 2), (ny * nid + 1):(ny * nid * 2)] / vf)
R.svd.h <- svd(R[(ny * nid * 2 + 1):(ny * nid * 3), (ny * nid * 2 + 1):(ny * nid * 3)] / vf)

R.root.u <- R.svd.u$u %*% diag(sqrt(R.svd.u$d))
R.root.v <- R.svd.v$u %*% diag(sqrt(R.svd.v$d))
R.root.h <- R.svd.h$u %*% diag(sqrt(R.svd.h$d))

epsilon.obs.meta <- matrix(NA, nrow = ny * nid, ncol = nobs)
for(i in 1:nobs){
  epsilon.obs.meta[, i] <- rnorm(ny * nid, 0, 1)
}

epsilon.obs.u <- R.root.u %*% epsilon.obs.meta
epsilon.obs.v <- R.root.v %*% epsilon.obs.meta
epsilon.obs.h <- R.root.h %*% epsilon.obs.meta

epsilon.obs.en.meta <- array(rnorm(ny * nid * n * nobs, 0, 1), c(ny * nid, n, nobs)) * sqrt(vf)
epsilon.obs.en <- array(NA, c(xobs, n, nobs))
for(i in 1:nobs){
  epsilon.obs.en[1:(ny * nid), , i] <- R.root.u %*% epsilon.obs.en.meta[, , i]
  epsilon.obs.en[(ny * nid + 1):(ny * nid * 2), , i] <- R.root.v %*% epsilon.obs.en.meta[, , i]
  epsilon.obs.en[(ny * nid * 2 + 1):(ny * nid * 3), , i] <- R.root.h %*% epsilon.obs.en.meta[, , i]
}
y.o[1:(ny * nid), ] <- so[1:(ny * nid), ] + epsilon.obs.u
y.o[(ny * nid + 1):(ny * nid * 2), ] <- so[(ny * nid + 1):(ny * nid * 2), ] + epsilon.obs.v
y.o[(ny * nid * 2 + 1):(ny * nid * 3), ] <- so[(ny * nid * 2 + 1):(ny * nid * 3), ] + epsilon.obs.h
for(i in 1:n) {
  y.f.en[, i, ] <- y.o + epsilon.obs.en[, i, ]
}
rm(R.svd.u, R.svd.v, R.svd.h,
   epsilon.obs.en, so, epsilon.obs.meta, epsilon.obs.en.meta, 
   epsilon.obs.u, epsilon.obs.v, epsilon.obs.h, R.root.u, R.root.v, R.root.h)

## Perturbed model state ####
x.en.bar = matrix(NA, ndim, S)
x.en.temp = matrix(NA, ndim, n)
epsilon.init = matrix(rnorm(ngrid * n, 0, sigma.i), nrow = ngrid)
H0 <- matrix(1, nrow = ny + 2, ncol = nx + 2)
U0 <- matrix(0, nrow = ny + 2, ncol = nx + 2)
V0 <- matrix(0, nrow = ny + 2, ncol = nx + 2)
for(j in 1:n){
  x.en.temp[(ngrid * 2 + 1):ndim, j] <- x.t[(ngrid * 2 + 1):ndim, 1] + epsilon.init[, j] # m ensemble
  H0[2:(ny + 1), 2:(nx + 1)] <- array(x.en.temp[(2 * ngrid + 1):(3 * ngrid), j], dim = c(ny, nx))
  H0[, 1] <- H0[, nx + 1]
  H0[, nx + 2] <- H0[, 2]
  H0[1, ] <- H0[2, ]
  H0[ny + 2, ] <- H0[ny + 1, ]
  k <- 2:(ny + 1)
  l <- 2:(nx + 1)
  U0[k, l] <- -g / f * (H0[k + 1, l] - H0[k - 1, l]) / (2 * dy)
  V0[k, l] <- g / f * (H0[k, l + 1] - H0[k, l - 1]) / (2 * dx)
  x.en.temp[1:(2 * ngrid), j] <- c(as.vector(U0[2:(ny + 1), 2:(nx + 1)]), as.vector(V0[2:(ny + 1), 2:(nx + 1)]))
}
x.en.bar[, 1] <- apply(x.en.temp, 1, mean)
rm(epsilon.init, H0, U0, V0)

## Output list ####
num = 0
cv = 0  ## convergence
vec.inflation.factor = rep(NA, nobs)
vec.objective.likelihood = rep(NA, nobs) 
vec.iteration.number = rep(NA, nobs)
vec.banding.bandwidth = matrix(NA, 6, nobs)
assimilate.errors = rep(NA, S)  ## RMSE at each step
assimilate.errors[1] = RMSE(x.t[, 1], x.en.bar[, 1])
for(t in 1:(S - 1)){
  ## One-step Forecast
  for(j in 1:n){
    x.en.temp[, j] <- SWE_Forecast(x.en.temp[, j], opt)
  }
  x.en.bar[, t + 1] <- apply(x.en.temp, 1, mean)
  
  ## Analyse if there exist observations
  if(t.o[t + 1] == 1){
    num <- num + 1
    
    x.f <- x.en.temp
    x.f.bar <- apply(x.f, 1, mean)
    d.f <- y.f.en[, , num] - H %*% x.f
    d.f.bar <- apply(d.f, 1, mean)
    
    # ### standard ####
    # Z = (x.f - x.f.bar) / sqrt(n - 1)
    # HZ = H %*% Z
    # K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
    # mulo = prod(1 + svd(R.inv.root %*% HZ)$d ^ 2)
    # obj.MLE = logdetR + log(mulo) + t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
    # 
    # inflation.fator = NA
    # objective.likelihood = obj.MLE
    # iteration.number = NA
    # banding.bandwidth = NA
    
    ### inflation & iteration ####
    Z = (x.f - x.f.bar) / sqrt(n - 1)
    HZ = H %*% Z
    D = svd(R.inv.root %*% HZ)$d
    log.likelihood = function(lambda){
      return(log(HPHR.det(lambda, D, R)) +
               t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
    }
    likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
    lambda0 = lambda = likelihood.optimize$minimum
    obj.MLE = likelihood.optimize$objective
    K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)

    # for(k in 1:10){
    #   x.a = x.f + K %*% d.f
    #   x.a.bar = apply(x.a, 1, mean)
    #   obj.MLE.old = obj.MLE
    # 
    #   Z = (x.f - x.a.bar) / sqrt(n - 1)
    #   HZ = H %*% Z
    #   D = svd(R.inv.root %*% HZ)$d
    #   log.likelihood = function(lambda){
    #     return(log(HPHR.det(lambda, D, R)) +
    #              t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
    #   }
    #   likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
    #   lambda = likelihood.optimize$minimum
    #   obj.MLE = likelihood.optimize$objective
    #   K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
    # 
    #   if(obj.MLE.old - obj.MLE < 1){
    #     break
    #   }
    # }

    x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar

    inflation.fator = lambda0
    objective.likelihood = obj.MLE
    iteration.number = NA
    banding.bandwidth = NA
    
    # ### inflation+taper ####
    # Z = (x.f - x.f.bar) / sqrt(n - 1)
    # HZ = H %*% Z
    # D = svd(R.inv.root %*% HZ)$d
    # log.likelihood = function(lambda){
    #   return(log(HPHR.det(lambda, D, R)) +
    #            t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
    # }
    # likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
    # lambda = likelihood.optimize$minimum
    # 
    # x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar  ## before all after ?
    # Z = sqrt(lambda) * Z
    # 
    # obj.taper = cal_taper_cov_SWE(t(Z), mat.dist, GC)
    # Z = obj.taper$Z
    # HZ = H %*% Z
    # D = svd(R.inv.root %*% HZ)$d
    # obj.MLE = log(HPHR.det(1, D, R)) +
    #   t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
    # K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
    # 
    # inflation.fator = lambda
    # objective.likelihood = obj.MLE
    # iteration.number = NA
    # banding.bandwidth = obj.taper$bandwidth
    
    ### Update ####
    x.en.temp = x.f + K %*% d.f
    x.en.bar[, t + 1] = apply(x.en.temp, 1, mean)
    vec.inflation.factor[num] = inflation.fator
    vec.objective.likelihood[num] = objective.likelihood
    vec.iteration.number[num] = iteration.number
    vec.banding.bandwidth[, num] = banding.bandwidth
    print(paste('Step', t, 'error', round(assimilate.errors[t], 3)))
  }
  assimilate.errors[t + 1] = RMSE(x.en.bar[, t + 1], x.t[, t + 1])
  
  if(assimilate.errors[t + 1] > 20 | is.na(assimilate.errors[t + 1])){
    cv = 1
    break
  }
}
t2 = Sys.time()
print(t2 - t1)
plot(assimilate.errors)
sqrt(mean(assimilate.errors[2941:5820] ^ 2))
sqrt(mean(assimilate.errors[5821:8700] ^ 2))
analyse.taper = x.en.bar
error.taper = assimilate.errors
save.image('SWE_refine_inflation.Rdata')
