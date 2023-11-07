rm(list = ls())
setwd("/disk/home/sunhaoxuan/EnKF")
library(tidyverse)
library(foreach)
library(doParallel)
library(Matrix)
source('Code/Functions/L96_Function_Oct23.R')

banding = function(z){
  if(z <= 1){
    return(1)
  }else{
    return(0)
  }
}

GC = function(z){
  z = z * 2
  if(z >= 0 & z <= 1){
    return(1 - 5 * z ^ 2 / 3 + 5 * z ^ 3 / 8 + z ^ 4 / 2 - z ^ 5 / 4)
  }else if(z > 1 & z <= 2){
    return(- 2 / 3 / z + 4 - 5 * z + 5 * z ^ 2 / 3 + 5 * z ^ 3 / 8 - z ^ 4 / 2 + z ^ 5 / 12)
  }else{
    return(0)
  }
}

cal_bandwidth_L96 = function(x, taper, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  coor = 1:p
  
  N_n = function(k){
    tmp = 0
    for(i in 1:(p - 1)){
      sigma_ii = sqrt(sum(x.scale[i, ] ^ 2))
      for(j in i:p){
        grid.distance = min(coor[j] - coor[i], p - coor[j] + coor[i])
        if(grid.distance <= k){
          sigma_ij = sum(x.scale[i, ] * x.scale[j, ])
          sigma_jj = sqrt(sum(x.scale[j, ] ^ 2))
          grid.weight = taper(grid.distance / k)
          tmp = tmp + 
            (grid.weight ^ 2 - 2 * grid.weight) * sigma_ij +
            (grid.weight ^ 2) * sigma_ii * sigma_jj /  n
        }
      }
    }
    return(tmp)
  }
  
  opt = optimise(N_n, c(0, p / 2), tol = 1) ## round ?
  bw = opt$minimum
  return(bw)
}

cal_mapping_L96_GC = function(x, bw, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  
  coor = 1:p
  
  obj.svd = svd(x.scale)
  u = obj.svd$u
  d = obj.svd$d
  d[which(d < 1e-5)] = 1e-5
  x.scale.inv = t(x.scale) %*% u %*% diag(1 / d ^ 2) %*% t(u)
  
  mapping.matrix = matrix(0, n, n)
  for(l in 1:p){
    Sigma2_ll = sum(x.scale[l, ] ^ 2)
    for(i in 1:n){
      mapping.matrix[i, i] = mapping.matrix[i, i] + x.scale.inv[i, l] ^ 2 * Sigma2_ll
    }
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
        mapping.matrix[i, j] = mapping.matrix[j, i] = mapping.matrix[i, j] + x.scale.inv[i, l] * x.scale.inv[j, l] * Sigma2_ll
      }
    }
  }
  for(l in 1:(p - 1)){
    for(m in (l + 1):p){
      grid.distance = min(coor[j] - coor[i], p - coor[j] + coor[i])
      if(grid.distance <= bw){
        Sigma2_lm = sum(x.scale[l, ] * x.scale[m, ]) * GC(grid.distance / k)
        for(i in 1:n){
          mapping.matrix[i, i] = mapping.matrix[i, i] + 2 * x.scale.inv[i, l] * x.scale.inv[i, m] * Sigma2_lm 
        }
        for(i in 1:(n - 1)){
          for(j in (i + 1):n){
            mapping.matrix[i, j] = mapping.matrix[j, i] = mapping.matrix[i, j] + (x.scale.inv[i, l] * x.scale.inv[j, m] + x.scale.inv[i, m] * x.scale.inv[j, l]) * Sigma2_lm
          }
        }
      }
    }
  }
  
  mapping.matrix.eigen = eigen(mapping.matrix)
  mapping.matrix.eigen.value = mapping.matrix.eigen$values
  mapping.matrix.eigen.vector = mapping.matrix.eigen$vectors
  mapping.matrix.eigen.value[which(mapping.matrix.eigen.value < 1e-5)] = 1e-5
  z = x.scale %*% mapping.matrix.eigen.vector %*% diag(sqrt(mapping.matrix.eigen.value))
  
  return(z)
}

dgTMat_L96_GC = function(x, bw, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  coor = 1:p
  
  banding.covariance = Matrix(0, p, p, sparse = T)
  for(i in 1:p){
    banding.covariance[i, i] = sum(x.scale[i, ] ^ 2)
  }
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      grid.distance = min(coor[j] - coor[i], p - coor[j] + coor[i])
      if(grid.distance <= bw){
        banding.covariance[i, j] = banding.covariance[j,  i] = sum(x.scale[i, ] * x.scale[j, ]) *  GC(grid.distance / bw)
      }
    }
  }
  
  return(banding.covariance)
}


f = 8
opt = list(
  p = 100,
  S = 2000,
  F.true = 8,
  h = 0.05,
  p.o = seq(1, 100, 1),
  S.o = seq(4, 2000, 4),
  sigma = 1,
  rho = 0.5,
  n = 30
)

state = L96_State(opt)

x.t = state$true
y.o = state$observation

p = opt$p
S = opt$S
F.true = opt$F.true
h = opt$h
sigma.tr = opt$sigma.tr
p.o = opt$p.o
S.o = opt$S.o
sigma = opt$sigma
rho = opt$rho
n = opt$n
F.set = f
sigma.in = f / 40

nsite = length(p.o)
nobs = length(S.o)
ti = seq(h, h * S, h)

H = matrix(0, nsite, p)
# H = Matrix(0, nsite, p, sparse = T)
for(k in 1:nsite){
  H[k, p.o[k]] = 1
}

R = matrix(0, nsite, nsite)
# R = Matrix(0, nsite, nsite, sparse = T)
for(i in 1:nsite){
  for(j in 1:nsite){
    R[i,j] = sigma ^ 2 * rho ^ min(abs(i - j), nsite - abs(i-j))
  }
}
R.inv = solve(R)
R.ev = eigen(R)
R.inv.root = R.ev$vectors %*% diag(R.ev$values ^ (-0.5)) %*% t(R.ev$vectors)
R.log.sum = sum(log(R.ev$values))

num = 0 ## observations
errs = rep(NA, S)  ## RMSE at each step
x.en.bar = matrix(NA, p, S)
x.en.temp = x.t[, 1] + t(MASS::mvrnorm(n, rep(0, p), diag(p) * sigma.in))
x.en.bar[, 1] = apply(x.en.temp, 1, mean)
errs[1] = RMSE(x.t[, 1], x.en.bar[, 1])
for(i in 1:(S - 1)){
  ## Update by RK4 
  for(j in 1:n){
    x.en.temp[, j] = RK4(x.en.temp[, j], ti[i], ti[i + 1] - ti[i], p, F.set)
  }
  x.en.bar[, i + 1] = apply(x.en.temp, 1, mean)
  
  ## Analyse if there exist observations
  if((i + 1) %in% S.o){
    num = num + 1
    
    x.f = x.en.temp
    x.f.bar = apply(x.f, 1, mean)
    d.f = y.o[p.o, num]  - H %*% x.f + t(MASS::mvrnorm(n, rep(0, nsite), R))
    d.f.bar = apply(d.f, 1, mean)
    
    # z = (x.f - x.f.bar) / sqrt(n - 1)
    # HZ = H %*% z
    # D = svd(R.inv.root %*% HZ)$d
    # log.likelihood = function(lambda){
    #   return(log(HPHR.det(lambda, D, R)) +
    #            t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
    # }
    # likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
    # lambda = likelihood.optimize$minimum
    # obj.MLE = likelihood.optimize$objective
    # K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
    # x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar

    bw = cal_bandwidth_L96_GC(x.f)
    P.hat = dgTMat_L96_GC(x.f, bw)
    P.eigen = eigen(P.hat)
    P.eigen.value = P.eigen$values
    P.eigen.vector = P.eigen$vectors
    z = t(P.eigen.vector[1:n, ]) %*% diag(P.eigen.value[1:n] ^ 0.5)
    HZ = H %*% z
    D = svd(R.inv.root %*% HZ)$d
    obj.MLE = log(HPHR.det(1, D, R)) +
      t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
    K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
    
    # bw = cal_bandwidth_L96_GC(x.f)
    # P.hat = dgTMat_L96_GC(x.f, bw)
    # HPHR.inv = solve(H %*% P.hat %*% t(H) + R + 1e-5 * diag(p))
    # K = P.hat %*% t(H) %*% HPHR.inv
    # obj.MLE = log(det(H %*% P.hat %*% t(H) + R)) + (t(d.f.bar) %*% HPHR.inv %*% d.f.bar)[1, 1]
    
    x.a = x.f + K %*% d.f
    x.a.bar = apply(x.a, 1, mean)
    for(k in 1:10){
      obj.MLE.old = obj.MLE
      
      # bw = cal_bandwidth_L96_GC(x.f)
      P.hat = dgTMat_L96_GC(x.f, bw, mu = x.a.bar)
      P.eigen = eigen(P.hat)
      P.eigen.value = P.eigen$values
      P.eigen.vector = P.eigen$vectors
      z = t(P.eigen.vector[1:n, ]) %*% diag(P.eigen.value[1:n] ^ 0.5)
      HZ = H %*% z
      D = svd(R.inv.root %*% HZ)$d
      obj.MLE = log(HPHR.det(1, D, R)) +
        t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
      K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
      
      if(obj.MLE.old - obj.MLE < 1){
        break
      }else{
        x.a = x.f + K %*% d.f
        x.a.bar = apply(x.a, 1, mean)
      }
    }

    x.en.temp = x.f + K %*% d.f
    x.en.bar[, i + 1] = apply(x.en.temp, 1, mean)
  }
  errs[i + 1] = RMSE(x.t[, i + 1], x.en.bar[, i + 1])
  print(paste(i, errs[i + 1]))
  
  ## Large error lead to divergence
  if(errs[i + 1] > 10 | is.nan(errs[i + 1])){
    cv = 1
    break
  }
}
# print(sqrt(mean(errs[1001:2000] ^ 2)))
