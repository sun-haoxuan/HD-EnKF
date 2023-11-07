## Test ####
rm(list = ls())
library(foreach)
library(Rcpp)

ngrid.x = 10
ngrid.y = 10
coor = data.frame(
  x = rep(1:ngrid.x, times = ngrid.y),
  y = rep(1:ngrid.y, each = ngrid.x)
)
p = ngrid.x * ngrid.y
n = 100
theta = 1 / 0.9
Sigma = diag(p)

## element-wise
for(i in 1:(p - 1)){
  for(j in i:p){
    grid.distance = sqrt(sum((coor[i, ] - coor[j, ]) ^ 2))
    Sigma[i, j] = Sigma[j, i] = theta ^ (-grid.distance)
  }
}
png('Eg_Sigma.png', width = 600, height = 600)
corrplot::corrplot(Sigma, method = 'color', is.corr = F)
dev.off()

corrplot::corrplot(Sigma, method = 'color', is.corr = F)
eg = eigen(Sigma)
A1 = eg$values
A1[which(A1 < 0)] = 0
Sigma.sqrt = eg$vectors %*% diag(A1 ^ 0.5) %*% t(eg$vectors)
x = Sigma.sqrt %*% matrix(rnorm(p * n, 0, 1), p, n)
x = (x - apply(x, 1, mean)) / sqrt(n - 1)

band.grids = function(bw, p, grids){
  band.grids = diag(1, p)
  for(i in 1:(p - 1)){
    for(j in i:p){
      coor1 = grids[i, ]
      coor2 = grids[j, ]
      x.dist = min(abs(coor1[1] - coor2[1]), ngrid.x -abs(coor1[1] - coor2[1]))
      y.dist = min(abs(coor1[2] - coor2[2]), ngrid.y -abs(coor1[2] - coor2[2]))
      grid.dist = sqrt(x.dist ^ 2 + y.dist ^ 2)
      if(grid.dist <= bw){
        band.grids[i, j] = band.grids[j, i] = 1
      }
    }
  }
  return(band.grids)
}

bw = 3
b.cal = band.grids(bw, p, coor)
png('Eg_Sigma_band.png', width = 600, height = 600)
corrplot::corrplot(Sigma * b.cal, method = 'color', is.corr = F)
dev.off()
png('Eg_Sigma_hat.png', width = 600, height = 600)
corrplot::corrplot(x %*% t(x) * b.cal, method = 'color', is.corr = F)
dev.off()

obj.svd = svd(x)
u = obj.svd$u
d = obj.svd$d
d[which(d < 1e-5)] = 1e-5
v = obj.svd$v
x.inv = t(x) %*% u %*% diag(1 / d ^ 2) %*%t(u)

t0 = Sys.time()
A = matrix(0, n, n)
for(k in 1:p){
  for(l in 1:p){
    if(b.cal[k, l] != 0){
      Sigma_kl = sum(x[k, ] * x[l, ])
      for(i in 1:n){
        for(j in 1:n){
          A[i, j] = A[i, j] + x.inv[i, k] * x.inv[j, l] * Sigma_kl
        }
      }
    }
  }
}
print(Sys.time() - t0)
sigma1 = x %*% A %*% t(x)

eg = eigen(A)
A1 = eg$vectors
A2 = eg$values
A2[which(A2 < 1e-5)] = 1e-5
z = x %*% A1 %*% diag(sqrt(A2))
png('Eg_Sigma_tilde.png', width = 600, height = 600)
corrplot::corrplot(z %*% t(z), method = 'color', is.corr = F)
dev.off()

cal_bandwidth = function(x, coor, lb, ub){
  p = nrow(x)

  M_n = function(k){
    tmp = 0
    for(i in 1:(p - 1)){
      sigma_ii = sqrt(mean(x[i, ] ^ 2))
      for(j in i:p){
        coor1 = coor[i, ]
        coor2 = coor[j, ]
        grid.distance = sqrt(sum((coor2 - coor1) ^ 2))
        if(grid.distance > k){
          tmp = tmp + mean(x[i, ] * x[j, ]) / p
        }else{
          tmp = tmp + sigma_ii * sqrt(mean(x[j, ] ^ 2)) / p / n
        }
      }
    }
    return(tmp)
  }
  
  opt = optimise(M_n, c(lb, ub), tol = 1)
  return(opt$minimum)
}

t1 = Sys.time()
cal_bandwidth(x, coor, 1, 10)
t2 = Sys.time()
print(t2 - t1)

cal_bandwidth = function(x, grids, n.split){
  p = nrow(x)
  n = ncol(x)
  
  set.seed(1234)
  n1 = n * (1 - 1 / log(n))
  index.list = foreach(s = 1:n.split) %do% {sample(1:n, n1)}
  
  # risk.F.norm = function(bw){
  #   tmp = 0
  #   for(i in 1:(p - 1)){
  #     for(j in i:p){
  #       # print(paste(i, j))
  #       coor1 = grids[i, ]
  #       coor2 = grids[j, ]
  #       x.dist = min(abs(coor1[1] - coor2[1]), 30 -abs(coor1[1] - coor2[1]))
  #       y.dist = min(abs(coor1[2] - coor2[2]), 30 -abs(coor1[2] - coor2[2]))
  #       grid.dist = sqrt(x.dist ^ 2 + y.dist ^ 2)
  # 
  #       sigma.ij.2 = foreach(s = 1:n.split, .combine = 'c') %do%
  #         {
  #           sum(x[i, -index.list[[s]]] * x[j, -index.list[[s]]])
  #         }
  #       if(grid.dist > bw){
  #         tmp = tmp + sum((sigma.ij.2) ^ 2)
  #       }else{
  #         sigma.ij.1 = foreach(s = 1:n.split, .combine = 'c') %do%
  #           {
  #             sum(x[i, index.list[[s]]] * x[j, index.list[[s]]])
  #           }
  #         tmp = tmp + sum((sigma.ij.1 - sigma.ij.2) ^ 2) / n.split / p
  #       }
  #     }
  #   }
  #   return(tmp)
  # }
  
  risk.F.norm = function(bw){
    grid.pair = list()
    for(i in 1:(p - 1)){
      for(j in i:p){
        coor1 = grids[i, ]
        coor2 = grids[j, ]
        x.dist = min(abs(coor1[1] - coor2[1]), ngrid.x -abs(coor1[1] - coor2[1]))
        y.dist = min(abs(coor1[2] - coor2[2]), ngrid.y -abs(coor1[2] - coor2[2]))
        grid.dist = sqrt(x.dist ^ 2 + y.dist ^ 2)
        if(grid.dist < bw){
          grid.pair[['i']] = c(grid.pair[['i']], i)
          grid.pair[['j']] = c(grid.pair[['j']], j)
        }
      }
    }
    grid.pair = as.data.frame(grid.pair)
    
    tmp = 0
    for(k in 1:nrow(grid.pair)){
      i = grid.pair$i[k]
      j = grid.pair$j[k]
      sigma.ij.1 = foreach(s = 1:n.split, .combine = 'c') %do%
        {
          sum(x[i, index.list[[s]]] * x[j, index.list[[s]]])
        }
      sigma.ij.2 = foreach(s = 1:n.split, .combine = 'c') %do%
        {
          sum(x[i, -index.list[[s]]] * x[j, -index.list[[s]]])
        }
      tmp = tmp + sum((sigma.ij.1 - sigma.ij.2) ^ 2) / n.split / p
    }
    return(tmp)
  }
  
  opt.risk = optimise(risk.F.norm, c(0, min(ngrid.x, ngrid.y) / 2), tol = 1)
  bw = opt.risk$minimum
  
  return(bw)
}




band.condition = function(i, j){
  coor1 = coor[i, ]
  coor2 = coor[j, ]
  x.dist = min(abs(coor1[1] - coor2[1]), ngrid.x -abs(coor1[1] - coor2[1]))
  y.dist = min(abs(coor1[2] - coor2[2]), ngrid.y -abs(coor1[2] - coor2[2]))
  grid.dist = sqrt(x.dist ^ 2 + y.dist ^ 2)
  return(grid.dist)
}

## L96 ####
rm(list = ls())
setwd("/disk/home/sunhaoxuan/EnKF")
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)

source('Code/Functions/band_Mat.R')
source('Code/Functions/cal_HPHR.R')
source('Code/Functions/cal_RMSE.R')
source('Code/Functions/cal_threshold.R')
source('Code/Functions/eig_adj.R')
source('Code/Functions/L96_State.R')
source('Code/Functions/L96_Analyse.R')

cal_bandwidth_L96 = function(x){
  p = nrow(x)
  n = ncol(x)
  x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  coor = matrix(c(cos(1:p * 2 * pi / p), sin(1:p * 2 * pi / p)), p, 2)
  
  M_n = function(k){
    tmp = 0
    for(i in 1:(p - 1)){
      sigma_ii = sqrt(sum(x.scale[i, ] ^ 2))
      for(j in i:p){
        coor1 = coor[i, ]
        coor2 = coor[j, ]
        grid.distance = sqrt(sum((coor2 - coor1) ^ 2))
        if(grid.distance > k){
          tmp = tmp + sum(x.scale[i, ] * x.scale[j, ]) / p
        }else{
          tmp = tmp + sigma_ii * sqrt(sum(x.scale[j, ] ^ 2)) / p / n
        }
      }
    }
    return(tmp)
  }
  
  opt = optimise(M_n, c(0, 2), tol = 1)
  return(opt$minimum)
}

cal_mapping_L96 = function(x, bw){
  p = nrow(x)
  n = ncol(x)
  x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  coor = matrix(c(cos(1:p * 2 * pi / p), sin(1:p * 2 * pi / p)), p, 2)
  
  obj.svd = svd(x.scale)
  u = obj.svd$u
  d = obj.svd$d
  d[which(d < 1e-5)] = 1e-5
  x.scale.inv = t(x.scale) %*% u %*% diag(1 / d ^ 2) %*% t(u)
  
  mapping.matrix = matrix(0, n, n)
  for(k in 1:p){
    for(l in 1:p){
      grid.distance = sqrt(sum((coor[k, ] - coor[l, ]) ^ 2))
        if(grid.distance <= bw){
          Sigma_kl = sum(x.scale[k, ] * x.scale[l, ])
          for(i in 1:n){
            for(j in 1:n){
              mapping.matrix[i, j] = mapping.matrix[i, j] + x.scale.inv[i, k] * x.scale.inv[j, l] * Sigma_kl
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

set.seed(1234)
opt = list(
  p = 40,
  S = 2000,
  F.true = 8,
  h = 0.05,
  sigma.tr = 0,
  p.o = seq(1, 40, 1),
  S.o = seq(4, 2000, 4),
  sigma = 1,
  rho = 0.5,
  n = 30,
  sigma.in = 0.1,
  F.set = 12,
  method = 'proposed',
  eig_adj = F,
  iter_update = T,
  is.circ = F
)
state = L96_State(opt)

### Initialize ####
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
F.set = opt$F.set
n = opt$n
sigma.in = opt$sigma.in
method = opt$method
eig_adj = opt$eig_adj
iter_update = opt$iter_update
is.circ = opt$circ

nsite = length(p.o)
nobs = length(S.o)
ti = seq(h, h * S, h)

### Observation matrix ####
H = matrix(0, nsite, p)
for(k in 1:nsite){
  H[k, p.o[k]] = 1
}

### Observation Error Covariance ####
R = matrix(0, nsite, nsite)
for(i in 1:nsite){
  for(j in 1:nsite){
    R[i,j] = sigma ^ 2 * rho ^ min(abs(i - j), nsite - abs(i-j))
  }
}
R.inv = solve(R)
R.ev = eigen(R)
R.inv.root = R.ev$vectors %*% diag(R.ev$values ^ (-0.5)) %*% t(R.ev$vectors)

### Output list ####
num = 0 ## observations
errs = rep(NA, S)  ## RMSE at each step
x.en.bar = matrix(NA, p, S)
x.en.temp = x.t[, 1] + mvrnorm(n, rep(0, p), diag(p) * sigma.in) %>% t()
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
    d.f = y.o[p.o, num]  - H %*% x.f
    d.en = d.f + t(mvrnorm(n, rep(0, nsite), R))
    d.f.bar = apply(d.en, 1, mean)
    
    # ## inflation
    # Z.f = 1 / sqrt(n - 1) * (x.f - x.f.bar)
    # HZ = H %*% Z.f
    # D = svd(R.inv.root %*% HZ)$d
    # ll = function(lambda){
    #   return(log(HPHR.det(lambda, D, R)) +
    #            t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
    # }
    # opt = optimize(ll, c(0.1, 10))
    # lambda = opt$minimum
    # K = lambda * Z.f %*% t(Z.f) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
    # x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
    
    ## banding
    bw = cal_bandwidth_L96(x.f)
    z = cal_mapping_L96(x.f, bw)
    HZ = H %*% z
    K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
    
    # ## banding + inflation
    # bw = cal_bandwidth_L96(x.f)
    # z = cal_mapping_L96(x.f, bw)
    # HZ = H %*% z
    # D = svd(R.inv.root %*% HZ)$d
    # ll = function(lambda){
    #   return(log(HPHR.det(lambda, D, R)) +
    #            t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
    # }
    # opt = optimize(ll, c(0.1, 10))
    # lambda = opt$minimum
    # K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
    # x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar

    x.en.temp = x.f + K %*% d.en
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
sqrt(mean(errs[1001:2000] ^ 2))
