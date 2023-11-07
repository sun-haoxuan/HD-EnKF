rm(list = ls())
gc()
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

cal_mapping_L96 = function(x, bw, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  
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

param = list()
for(f in 4:12){
  for(b in 1:2){
    for(inflate in c(T, F)){
      for(iterative.update in c(T, F)){
        param[['F.set']] = c(param[['F.set']], f)
        param[['boot']] = c(param[['boot']], b)
        param[['inflate']] = c(param[['inflate']], inflate)
        param[['iterative.update']] = c(param[['iterative.update']], iterative.update)
      }
    }
  }
}
param = as.data.frame(param)


cores = 36
cl = makeCluster(cores)
registerDoParallel(cl, cores = cores)
result = foreach(i = 1:nrow(param), .combine = 'rbind') %dopar%
  {
    f = param$F.set[i]
    b = param$boot[i]
    inflate = param$inflate[i]
    iterative.update = param$iterative.update[i]
    
    set.seed(1234 + b)
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
      F.set = f,
      inflate = inflate,
      iterative.update = iterative.update
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
    iterative.update = opt$iterative.update
    inflate = opt$inflate
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
    R.log.sum = sum(log(R.ev$values))
    
    ### Output list ####
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
        d.f = y.o[p.o, num]  - H %*% x.f
        d.en = d.f + t(MASS::mvrnorm(n, rep(0, nsite), R))
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
        
        # ## banding
        # bw = cal_bandwidth_L96(x.f)
        # z = cal_mapping_L96(x.f, bw)
        # HZ = H %*% z
        # K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        
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
        
        ## proposed
        bw = cal_bandwidth_L96(x.f)
        z = cal_mapping_L96(x.f, bw)
        HZ = H %*% z
        
        if(inflate){
          D = svd(R.inv.root %*% HZ)$d
          ll = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          opt = optimize(ll, c(0.1, 10))
          lambda = opt$minimum
          x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
          K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
        }else{
          K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        }
        
        if(iterative.update){
          mulo = prod(1 + svd(R.inv.root %*% HZ)$d ^ 2)
          obj.MLE = R.log.sum + log(mulo) + t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
          
          x.a = x.f + K %*% d.en
          x.a.bar = apply(x.a, 1, mean)
          
          for(k in 1:10){
            obj.MLE.old = obj.MLE
            K.old = K
            
            z = cal_mapping_L96(x.f, bw, mu = x.a.bar)
            HZ = H %*% z
            if(inflate){
              D = svd(R.inv.root %*% HZ)$d
              ll = function(lambda){
                return(log(HPHR.det(lambda, D, R)) +
                         t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
              }
              opt = optimize(ll, c(0.1, 10))
              lambda = opt$minimum
              x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
              K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
            }else{
              K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
            }
            mulo = prod(1 + svd(R.inv.root %*% HZ)$d ^ 2)
            obj.MLE = R.log.sum + log(mulo) + t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }else{
              x.a = x.f + K %*% d.en
              x.a.bar = apply(x.a, 1, mean)
            }
          }
          obj.MLE = obj.MLE.old
          K = K.old
        }
        
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
    
    data.frame(
      f = f,
      boot = b,
      inflate = inflate,
      iterative.update = iterative.update,
      rmse = sqrt(mean(errs[1001:2000] ^ 2))
    )
  }
stopCluster(cl)
stopImplicitCluster()

result %>% 
  mutate(group = paste0('infl', inflate, '_iter', iterative.update)) %>% 
  group_by(group, f) %>% 
  summarise(rmse = min(rmse, na.rm = T)) %>% 
  ggplot() + geom_line(aes(x = f, y = rmse, color = group))