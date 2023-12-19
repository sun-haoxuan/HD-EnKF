L96_Analyse = function(state, option){
  
  require(MASS)
  
  ## Initialize ####
  x.t = state$true
  y.o = state$observation
  
  p = option$p
  S = option$S
  F.true = option$F.true
  h = option$h
  p.o = option$p.o
  S.o = option$S.o
  sigma = option$sigma
  rho = option$rho
  F.set = option$F.set
  n = option$n
  sigma.in = option$sigma.in
  method = option$method
  is.iter = option$is.iter
  interval = option$interval
  
  nsite = length(p.o)
  nobs = length(S.o)
  ti = seq(h, h * S, h)
  
  ## Observation matrix ####
  H = matrix(0, nsite, p)
  for(k in 1:nsite){
    H[k, p.o[k]] = 1
  }
  
  ## Observation Error Covariance ####
  R = matrix(0, nsite, nsite)
  for(i in 1:nsite){
    for(j in 1:nsite){
      R[i,j] = sigma ^ 2 * rho ^ min(abs(i - j), nsite - abs(i-j))
    }
  }
  R.inv = solve(R)
  R.ev = eigen(R)
  R.inv.root = R.ev$vectors %*% diag(R.ev$values ^ (-0.5)) %*% t(R.ev$vectors)
  
  ## distance matrix
  mat.dist = diag(0, p)
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      mat.dist[i, j] = mat.dist[j, i] = min(j - i, p - j + i)
    }
  }
  
  ## Output list ####
  cv = 0  ## convergence
  num = 0 ## observations
  k = 0 ## iteration
  assimilate.errors = rep(NA, S)  ## RMSE at each step
  x.en.bar = matrix(NA, p, S)
  x.en.temp = x.t[, 1] + t(MASS::mvrnorm(n, rep(0, p), diag(p) * sigma.in))
  x.en.bar[, 1] = apply(x.en.temp, 1, mean)
  
  ## Analyse State ####
  assimilate.errors[1] = RMSE(x.t[, 1], x.en.bar[, 1])
  vec.inflation.factor = rep(NA, nobs)
  vec.objective.likelihood = rep(NA, nobs) 
  vec.iteration.number = rep(NA, nobs)
  vec.banding.bandwidth = rep(NA, nobs)
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
      
      ### standard ####
      if(method == 'standard'){
        z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = NA
      }
      
      ### Oracle ####
      if(method == 'oracle'){
        z = (x.f - x.t[, i + 1]) / sqrt(n - 1)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = NA
      }
      
      ### Inflation ####
      if(method == 'inflation'){
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
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
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
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = NA
      }
      
      ### banding-inflation ####
      if(method == 'banding-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        k = NA
        if(is.iter){
          for(k in 1:10){
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
            obj.MLE.old = obj.MLE
            K.old = K
            
            Z = (x.f - x.a.bar) / sqrt(n - 1)
            Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
            HZ = H %*% Z
            D = svd(R.inv.root %*% HZ)$d
            log.likelihood = function(lambda){
              return(log(HPHR.det(lambda, D, R)) +
                       t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
            }
            likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
            lambda0 = lambda = likelihood.optimize$minimum
            obj.MLE = likelihood.optimize$objective
            K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### tapering-inflation ####
      if(method == 'tapering-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        k = NA
        if(is.iter){
          for(k in 1:10){
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
            obj.MLE.old = obj.MLE
            K.old = K
            
            Z = (x.f - x.a.bar) / sqrt(n - 1)
            Z = cal_taper_Z(t(Z), bw, mat.dist, tapering)
            HZ = H %*% Z
            D = svd(R.inv.root %*% HZ)$d
            log.likelihood = function(lambda){
              return(log(HPHR.det(lambda, D, R)) +
                       t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
            }
            likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
            lambda0 = lambda = likelihood.optimize$minimum
            obj.MLE = likelihood.optimize$objective
            K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### GC-inflation ####
      if(method == 'GC-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        k = NA
        if(is.iter){
          for(k in 1:10){
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
            obj.MLE.old = obj.MLE
            K.old = K
            
            Z = (x.f - x.a.bar) / sqrt(n - 1)
            Z = cal_taper_Z(t(Z), bw, mat.dist, GC)
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
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### banding-ev-p ####
      if(method == 'banding-ev-p'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        k = NA
        if(is.iter){
          for(k in 1:10){
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
            obj.MLE.old = obj.MLE
            K.old = K
            
            Z = (x.f - x.a.bar) / sqrt(n - 1)
            Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
            HZ = H %*% Z
            D = svd(R.inv.root %*% HZ)$d
            log.likelihood = function(lambda){
              return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                       t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
            }
            likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
            lambda = likelihood.optimize$minimum
            obj.MLE = likelihood.optimize$objective
            K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### tapering-ev-p ####
      if(method == 'tapering-ev-p'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        k = NA
        if(is.iter){
          for(k in 1:10){
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
            obj.MLE.old = obj.MLE
            K.old = K
            
            Z = (x.f - x.a.bar) / sqrt(n - 1)
            Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = p)
            HZ = H %*% Z
            D = svd(R.inv.root %*% HZ)$d
            log.likelihood = function(lambda){
              return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                       t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
            }
            likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
            lambda0 = lambda = likelihood.optimize$minimum
            obj.MLE = likelihood.optimize$objective
            K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### GC-ev-p ####
      if(method == 'GC-ev-p'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective

        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        k = NA
        if(is.iter){
          for(k in 1:10){
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
            obj.MLE.old = obj.MLE
            K.old = K
            
            Z = (x.f - x.a.bar) / sqrt(n - 1)
            Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
            HZ = H %*% Z
            D = svd(R.inv.root %*% HZ)$d
            log.likelihood = function(lambda){
              return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                       t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
            }
            likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
            lambda0 = lambda = likelihood.optimize$minimum
            obj.MLE = likelihood.optimize$objective
            K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### Update ####
      x.en.temp = x.f + K %*% d.f
      x.en.bar[, i + 1] = apply(x.en.temp, 1, mean)
      vec.inflation.factor[num] = inflation.fator
      vec.objective.likelihood[num] = objective.likelihood
      vec.iteration.number[num] = iteration.number
      vec.banding.bandwidth[num] = banding.bandwidth
    }
    assimilate.errors[i + 1] = RMSE(x.t[, i + 1], x.en.bar[, i + 1])
    print(paste(i, assimilate.errors[i + 1]))
    
    ## Large error lead to divergence
    if(assimilate.errors[i + 1] > 10 | is.nan(assimilate.errors[i + 1])){
      cv = 1
      break
    }
  }
  print(sqrt(mean(assimilate.errors[1001:2000] ^ 2)))
  
  output = list(
    analyse = x.en.bar, 
    converge = cv, 
    error = assimilate.errors, 
    option = option,
    inflation = vec.inflation.factor,
    likelihood = vec.objective.likelihood,
    iteration = vec.iteration.number,
    bandwidth = vec.banding.bandwidth
  )
  
  return(output)
}