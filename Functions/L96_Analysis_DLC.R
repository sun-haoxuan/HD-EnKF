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
      
      ### Oracle-inflation ####
      if(method == 'oracle-inflation'){
        Z = (x.f - x.t[, i + 1]) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.t[, i + 1]) + x.t[, i + 1]
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = NA
      }
      
      ### Oracle-banding ####
      if(method == 'oracle-banding'){
        Z = (x.f - x.t[, i + 1]) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det2(1, D, R, sd(D))) + length(D) * log(sd(D))+
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
        x.f = (x.f - x.t[, i + 1]) + x.t[, i + 1]
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = bw
      }
      
      ### Oracle-banding-inflation ####
      if(method == 'oracle-banding-inflation'){
        Z = (x.f - x.t[, i + 1]) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.t[, i + 1]) + x.t[, i + 1]
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = bw
      }
      
      ### Oracle-GC ####
      if(method == 'oracle-GC'){
        Z = (x.f - x.t[, i + 1]) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det2(1, D, R, sd(D))) + length(D) * log(sd(D))+
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
        x.f = (x.f - x.t[, i + 1]) + x.t[, i + 1]
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = bw
      }
      
      ### Oracle-GC-inflation ####
      if(method == 'oracle-GC-inflation'){
        Z = (x.f - x.t[, i + 1]) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        x.f = sqrt(lambda0) * (x.f - x.t[, i + 1]) + x.t[, i + 1]
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = bw
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
      
      ### banding ####
      if(method == 'banding'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, min(n, p)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, min(n, p))
        
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
            obj.MLE = log(HPHR.det(1, D, R)) +
              t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, min(n, p)) %*% d.f.bar
            K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, min(n, p))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### tapering ####
      if(method == 'tapering'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, min(n, p)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, min(n, p))
        
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
            obj.MLE = log(HPHR.det(1, D, R)) +
              t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, min(n, p)) %*% d.f.bar
            K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, min(n, p))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### GC ####
      if(method == 'GC'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, min(n, p)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, min(n, p))
        
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
            obj.MLE = log(HPHR.det(1, D, R)) +
              t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, min(n, p)) %*% d.f.bar
            K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, min(n, p))
            
            if(obj.MLE.old - obj.MLE < 1){
              break
            }
          }
          obj.MLE = obj.MLE
          K = K.old
        }
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
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
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
        
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
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
        
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
      
      ### GC-inflation ####
      if(method == 'GC-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
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
      
      ### tapering-ev-p ####
      if(method == 'tapering-ev-p'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = p)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = p)
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
      
      ### GC-ev-p ####
      if(method == 'GC-ev-p'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
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
      
      ### banding-smoothBW ####
      if(method == 'banding-smoothBW'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)])
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      
      ### banding-smoothBW-1 ####
      if(method == 'banding-smoothBW-1'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) - 1
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      ### banding-smoothBW+1 ####
      if(method == 'banding-smoothBW+1'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) + 1
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      
      ### banding-smoothBW-2 ####
      if(method == 'banding-smoothBW-2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) - 2
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      ### banding-smoothBW+2 ####
      if(method == 'banding-smoothBW+2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) + 2
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      
      ### banding-smoothBW_split2 ####
      if(method == 'banding-smoothBW_split2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) / 2
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      ### banding-smoothBW_multi2 ####
      if(method == 'banding-smoothBW_multi2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) * 2
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      
      ### GC-smoothBW ####
      if(method == 'GC-smoothBW'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)])
        }
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
      
      ### GC-smoothBW-1 ####
      if(method == 'GC-smoothBW-1'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) - 1
        }
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
      ### GC-smoothBW+1 ####
      if(method == 'GC-smoothBW+1'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) + 1
        }
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
      
      ### GC-smoothBW-2 ####
      if(method == 'GC-smoothBW-2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) - 2
        }
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
      ### GC-smoothBW+2 ####
      if(method == 'GC-smoothBW+2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) + 2
        }
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
      
      ### GC-smoothBW_split2 ####
      if(method == 'GC-smoothBW_split2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) /2
        }
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
      ### GC-smoothBW_multi2 ####
      if(method == 'GC-smoothBW_multi2'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)]) * 2
        }
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
      
      ### banding-smoothBW ####
      if(method == 'banding-smoothBW'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)])
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding)
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
      
      ### tapering-smoothBW ####
      if(method == 'tapering-smoothBW'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)])
        }
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
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
        
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
      
      ### GC-smoothBW ####
      if(method == 'GC-smoothBW'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)])
        }
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
      
      ### banding-ev30 ####
      if(method == 'banding-ev30'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = 30)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = 30)
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
      
      ### tapering-ev30 ####
      if(method == 'tapering-ev30'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = 30)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = 30)
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
      
      ### GC-ev30 ####
      if(method == 'GC-ev30'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = 30)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = 30)
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
      
      ### banding-ev60 ####
      if(method == 'banding-ev60'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = 60)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = 60)
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
      
      ### tapering-ev60 ####
      if(method == 'tapering-ev60'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = 60)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = 60)
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
      
      ### GC-ev60 ####
      if(method == 'GC-ev60'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = 60)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = 60)
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
      
      ### banding-ev90 ####
      if(method == 'banding-ev90'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = 90)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = 90)
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
      
      ### tapering-ev90 ####
      if(method == 'tapering-ev90'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = 90)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = 90)
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
      
      ### GC-ev90 ####
      if(method == 'GC-ev90'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = 90)
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
            Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = 90)
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
      
      ### banding-approx ####
      if(method == 'banding-approx'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z_approx(t(Z), bw, mat.dist, banding)
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
            Z = cal_taper_Z_approx(t(Z), bw, mat.dist, banding)
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
      
      ### tapering-approx ####
      if(method == 'tapering-approx'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z_approx(t(Z), bw, mat.dist, tapering)
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
            Z = cal_taper_Z_approx(t(Z), bw, mat.dist, tapering)
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
      
      ### GC-approx ####
      if(method == 'GC-approx'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z_approx(t(Z), bw, mat.dist, GC)
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
            Z = cal_taper_Z_approx(t(Z), bw, mat.dist, GC)
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
      
      ### banding-only ####
      if(method == 'banding-only'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### tapering-only ####
      if(method == 'tapering-only'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, tapering, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, tapering, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### GC-only ####
      if(method == 'GC-only'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### GC-constant ####
      if(method == 'GC-constant'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        if(i < S / 2){
          bw = cal_taper_bandwidth(t(Z), mat.dist, GC, interval)
        }else{
          bw = mean(vec.banding.bandwidth[1:(nobs / 2)], trim = .1)
        }
        Z = cal_taper_Z(t(Z), bw, mat.dist, GC, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
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
              return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
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
      
      ### banding-inflation(multiplicative) ####
      if(method == 'banding-inflation(multiplicative)'){
        lambda0 = obj.MLE = k = bw = NA
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
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### banding-inflation(additive) ####
      if(method == 'banding-inflation(additive)'){
        lambda0 = obj.MLE = k = bw = NA
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, banding, interval)
        Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det2(1, D, R, sd(D))) + p * log(sd(D))+
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z(t(Z), bw, mat.dist, banding, ndim = p)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          obj.MLE = log(HPHR.det2(1, D, R, sd(D))) + p * log(sd(D))+
            t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### banding-inflation(both) ####
      if(method == 'banding-inflation(both)'){
        lambda0 = obj.MLE = k = bw = NA
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
          lambda0 = lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### banding-block ####
      if(method == 'banding-block'){
        lambda0 = obj.MLE = k = bw = NA
        n.b = n.b
        ol.b = ol.b
        p.b = p / n.b
        index.list = list()
        for(b in 1:n.b){
          index.tmp = (1 - ol.b):(p.b + ol.b) + p.b * (b - 1)
          index.tmp[which(index.tmp <= 0)] = index.tmp[which(index.tmp <= 0)] + p
          index.tmp[which(index.tmp > p)] = index.tmp[which(index.tmp > p)] - p
          weight.tmp = c((1:ol.b) / (ol.b + 1), rep(1, p.b), (ol.b:1) / (ol.b + 1))
          index.list[[b]] = data.frame(index = index.tmp, weight = weight.tmp)
        }
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        Z1 = matrix(0, p, p)
        bandwidth.vec = rep(NA, n.b)
        for(b in 1:n.b){
          index.b = index.list[[b]]$index
          weight.b = index.list[[b]]$weight
          Z.b = Z[index.b, ]
          bw.b = cal_taper_bandwidth(t(Z.b), mat.dist[index.b, index.b], banding, interval)
          Z.b = cal_taper_Z(t(Z.b), bw.b, mat.dist[index.b, index.b], banding, ndim = p.b + ol.b * 2)
          Z1[index.b, index.b] = Z1[index.b, index.b] + Z.b %*% diag(weight.b)
          bandwidth.vec[b] = bw.b
        }
        Z = Z1
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z1 = matrix(0, p, p)
          for(b in 1:n.b){
            index.b = index.list[[b]]$index
            weight.b = index.list[[b]]$weight
            Z.b = Z[index.b, ]
            bw.b = bandwidth.vec[b]
            Z.b = cal_taper_Z(t(Z.b), bw.b, mat.dist[index.b, index.b], banding, ndim = p.b + ol.b * 2)
            Z1[index.b, index.b] = Z1[index.b, index.b] + Z.b %*% diag(weight.b)
          }
          Z = Z1
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda0 = lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
          K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### tapering-block ####
      if(method == 'tapering-block'){
        lambda0 = obj.MLE = k = bw = NA
        n.b = n.b
        ol.b = ol.b
        p.b = p / n.b
        index.list = list()
        for(b in 1:n.b){
          index.tmp = (1 - ol.b):(p.b + ol.b) + p.b * (b - 1)
          index.tmp[which(index.tmp <= 0)] = index.tmp[which(index.tmp <= 0)] + p
          index.tmp[which(index.tmp > p)] = index.tmp[which(index.tmp > p)] - p
          weight.tmp = c((1:ol.b) / (ol.b + 1), rep(1, p.b), (ol.b:1) / (ol.b + 1))
          index.list[[b]] = data.frame(index = index.tmp, weight = weight.tmp)
        }
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        Z1 = matrix(0, p, p)
        bandwidth.vec = rep(NA, n.b)
        for(b in 1:n.b){
          index.b = index.list[[b]]$index
          weight.b = index.list[[b]]$weight
          Z.b = Z[index.b, ]
          bw.b = cal_taper_bandwidth(t(Z.b), mat.dist[index.b, index.b], tapering, interval)
          Z.b = cal_taper_Z(t(Z.b), bw.b, mat.dist[index.b, index.b], tapering, ndim = p.b + ol.b * 2)
          Z1[index.b, index.b] = Z1[index.b, index.b] + Z.b %*% diag(weight.b)
          bandwidth.vec[b] = bw.b
        }
        Z = Z1
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z1 = matrix(0, p, p)
          for(b in 1:n.b){
            index.b = index.list[[b]]$index
            weight.b = index.list[[b]]$weight
            Z.b = Z[index.b, ]
            bw.b = bandwidth.vec[b]
            Z.b = cal_taper_Z(t(Z.b), bw.b, mat.dist[index.b, index.b], tapering, ndim = p.b + ol.b * 2)
            Z1[index.b, index.b] = Z1[index.b, index.b] + Z.b %*% diag(weight.b)
          }
          Z = Z1
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda0 = lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
          K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### GC-block ####
      if(method == 'GC-block'){
        lambda0 = obj.MLE = k = bw = NA
        n.b = n.b
        ol.b = ol.b
        p.b = p / n.b
        index.list = list()
        for(b in 1:n.b){
          index.tmp = (1 - ol.b):(p.b + ol.b) + p.b * (b - 1)
          index.tmp[which(index.tmp <= 0)] = index.tmp[which(index.tmp <= 0)] + p
          index.tmp[which(index.tmp > p)] = index.tmp[which(index.tmp > p)] - p
          weight.tmp = c((1:ol.b) / (ol.b + 1), rep(1, p.b), (ol.b:1) / (ol.b + 1))
          index.list[[b]] = data.frame(index = index.tmp, weight = weight.tmp)
        }
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        Z1 = matrix(0, p, p)
        bandwidth.vec = rep(NA, n.b)
        for(b in 1:n.b){
          index.b = index.list[[b]]$index
          weight.b = index.list[[b]]$weight
          Z.b = Z[index.b, ]
          bw.b = cal_taper_bandwidth(t(Z.b), mat.dist[index.b, index.b], GC, interval)
          Z.b = cal_taper_Z(t(Z.b), bw.b, mat.dist[index.b, index.b], GC, ndim = p.b + ol.b * 2)
          Z1[index.b, index.b] = Z1[index.b, index.b] + Z.b %*% diag(weight.b)
          bandwidth.vec[b] = bw.b
        }
        Z = Z1
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z1 = matrix(0, p, p)
          for(b in 1:n.b){
            index.b = index.list[[b]]$index
            weight.b = index.list[[b]]$weight
            Z.b = Z[index.b, ]
            bw.b = bandwidth.vec[b]
            Z.b = cal_taper_Z(t(Z.b), bw.b, mat.dist[index.b, index.b], GC, ndim = p.b + ol.b * 2)
            Z1[index.b, index.b] = Z1[index.b, index.b] + Z.b %*% diag(weight.b)
          }
          Z = Z1
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det2(lambda, D, R, sd(D))) + p * log(sd(D))+
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda0 = lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
          K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### sum-eigen ####
      if(method == 'sum-eigen'){
        ### assimilation ####
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth(t(Z), mat.dist, taper, interval)
        Z = cal_taper_Z2(t(Z), bw, mat.dist, taper, sumeigen)
        vec.number.eigen[num] = ncol(Z)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        
        K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        k = NA
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z2(t(Z), bw, mat.dist, taper, 0.9)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
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