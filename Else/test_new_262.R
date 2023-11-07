## Test ####
rm(list = ls())
library(foreach)
library(Rcpp)

ngrid.x = 20
ngrid.y = 20
grids = data.frame(
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
    print(paste(i, j))
    coor1 = grids[i, ]
    coor2 = grids[j, ]
    x.dist = min(abs(coor1[1] - coor2[1]), ngrid.x -abs(coor1[1] - coor2[1]))
    y.dist = min(abs(coor1[2] - coor2[2]), ngrid.y -abs(coor1[2] - coor2[2]))
    Sigma[i, j] = Sigma[j, i] = theta ^ (-sqrt(x.dist ^ 2 + y.dist ^ 2))
  }
}

## row-wise
png('Eg_Sigma.png', width = 600, height = 600)
corrplot::corrplot(Sigma, method = 'color', is.corr = F)
dev.off()

eg = eigen(Sigma)
A1 = eg$values
A1[which(A1 < 0)] = 0
Sigma.sqrt = eg$vectors %*% diag(A1 ^ 0.5) %*% t(eg$vectors)
x = Sigma.sqrt %*% matrix(rnorm(p * n, 0, 1), p, n)

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

bw = 2
b.cal = band.grids(bw, p, grids)
corrplot::corrplot(b.cal, method = 'color', is.corr = F)
sigma = x %*% t(x) * b.cal
corrplot::corrplot(sigma, method = 'color', is.corr = F)


band.condition = function(i, j){
  coor1 = grids[i, ]
  coor2 = grids[j, ]
  x.dist = min(abs(coor1[1] - coor2[1]), ngrid.x -abs(coor1[1] - coor2[1]))
  y.dist = min(abs(coor1[2] - coor2[2]), ngrid.y -abs(coor1[2] - coor2[2]))
  grid.dist = sqrt(x.dist ^ 2 + y.dist ^ 2)
  return(grid.dist)
}

obj.svd = svd(x)
u = obj.svd$u
d = obj.svd$d
v = obj.svd$v
d[which(d < 1e-5)] = 1e-5
x.inv = t(x) %*% u %*% diag(1 / d ^ 2) %*%t(u)

t0 = Sys.time()
A0 = x.inv %*% sigma %*% t(x.inv)
A = matrix(0, n, n)
for(l in 1:p){
  for(m in 1:p){
    if(b.cal[l, m] != 0){
    # if(band.condition(k, l) <= bw){
      Sigma_lm = sum(x[l, ] * x[m, ])
      # for(i in 1:n){
      #   for(j in i:n){
      #     A[i, j] = A[j, i] = A[i, j] + x.inv[i, l] * x.inv[j, m] * Sigma_lm
      #   }
      # }
      for(i in 1:n){
        A[i, i] = A[i, i] + x.inv[i, l] * x.inv[i, m] * Sigma_lm
      }
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          A[i, j] = A[j, i] = A[i, j] + x.inv[i, l] * x.inv[j, m] * Sigma_lm
        }
      }
    }
  }
}
A = matrix(0, n, n)
for(l in 1:p){
  Sigma_ll = sum(x[l, ] ^ 2)
  for(i in 1:n){
    A[i, i] = A[i, i] + x.inv[i, l] ^ 2 * Sigma_ll
  }
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      A[i, j] = A[j, i] = A[i, j] + x.inv[i, l] * x.inv[j, l] * Sigma_ll
    }
  }
}
for(l in 1:(p - 1)){
  for(m in (l + 1):p){
    if(b.cal[l, m] != 0){
      Sigma_lm = sum(x[l, ] * x[m, ])
      for(i in 1:n){
        A[i, i] = A[i, i] + (x.inv[i, l] * x.inv[i, m] + x.inv[i, m] * x.inv[i, l]) * Sigma_lm
      }
      for(i in 1:(n - 1)){
        for(j in (i + 1):n){
          A[i, j] = A[j, i] = A[i, j] + (x.inv[i, l] * x.inv[j, m] + x.inv[i, m] * x.inv[j, l]) * Sigma_lm
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
sigma2 = z %*% t(z)
corrplot::corrplot(sigma2, method = 'color', is.corr = F)

## L96 ####
rm(list = ls())
setwd("E:/Project/EnKF")
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)

cal_bandwidth_L96 = function(x, p, n.split){
  p = nrow(x)
  n = ncol(x)
  
  grids = 1:p
  
  set.seed(1234)
  n1 = n * (1 - 1 / log(n))
  index.list = foreach(s = 1:n.split) %do% {sample(1:n, n1)}
  
  risk.F.norm = function(bw){
    grid.pair = list()
    for(i in 1:(p - 1)){
      for(j in i:p){
        coor1 = grids[i]
        coor2 = grids[j]
        grid.dist = min(abs(coor1 - coor2), p -abs(coor1 - coor2))
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
  
  opt.risk = optimise(risk.F.norm, c(0, p / 2), tol = 1)
  bw = opt.risk$minimum
  
  return(bw)
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
    
    x1 = x.f - apply(x.f, 1, mean)
    
    
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