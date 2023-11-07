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
  
  nsite = length(p.o)
  nobs = length(S.o)
  ti = seq(h, h * S, h)
  
  ## Observation matrix ####
  H = matrix(0, nsite, p)
  # H = Matrix(0, nsite, p, sparse = T)
  for(k in 1:nsite){
    H[k, p.o[k]] = 1
  }
  
  ## Observation Error Covariance ####
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
  
  ## Output list ####
  cv = 0  ## convergence
  num = 0 ## observations
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
      
      ### Wang & Bishop ####
      if(method == 'Wang & Bishop'){
        z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        lambda = (t(d.f.bar) %*% t(R.inv.root) %*% R.inv.root %*% d.f.bar - p)[1, 1] /
          sum(diag(R.inv.root %*% HZ %*% t(HZ) %*% t(R.inv.root)))
        lambda = max(lambda, 0.1)
        K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        obj.MLE = log(HPHR.det(lambda, D, R)) +
          t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = NA
      }
      
      ### Inflation ####
      if(method == 'inflation'){
        z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = NA
      }
      
      ### Inflation & Iteration ####
      if(method == 'inflation & iteration'){
        z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        x.a = x.f + K %*% d.f
        x.a.bar = apply(x.a, 1, mean)
        for(k in 1:10){
          obj.MLE.old = obj.MLE
          
          z = (x.f - x.a.bar) / sqrt(n - 1)
          HZ = H %*% z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }else{
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
          }
        }
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = NA
      }
      
      ### Banding ####
      if(method == 'banding'){
        bw = cal_bandwidth_L96(x.f)
        z = cal_mapping_L96(x.f, bw)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = bw
      }
      
      ### Banding & Inflation ####
      if(method == 'banding & inflation'){
        bw = cal_bandwidth_L96(x.f)
        z = cal_mapping_L96(x.f, bw)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = bw
      }
      
      ### Banding & Iteration ####
      if(method == 'banding & iteration'){
        bw = cal_bandwidth_L96(x.f)
        z = cal_mapping_L96(x.f, bw)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = z %*% t(z) %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.a = x.f + K %*% d.f
        x.a.bar = apply(x.a, 1, mean)
        for(k in 1:10){
          obj.MLE.old = obj.MLE
          
          z = cal_mapping_L96(x.f, bw, mu = x.a.bar)
          HZ = H %*% z
          D = svd(R.inv.root %*% HZ)$d
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
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### Banding & Inflation & Iteration ####
      if(method == 'banding & inflation & iteration'){
        bw = cal_bandwidth_L96(x.f)
        z = cal_mapping_L96(x.f, bw)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        x.a = x.f + K %*% d.f
        x.a.bar = apply(x.a, 1, mean)
        for(k in 1:10){
          obj.MLE.old = obj.MLE
          
          z = cal_mapping_L96(x.f, bw, mu = x.a.bar)
          HZ = H %*% z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }else{
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
          }
        }
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### inflation iteration then banding ####
      if(method == 'inflation iteration then banding'){
        z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        x.a = x.f + K %*% d.f
        x.a.bar = apply(x.a, 1, mean)
        for(k in 1:10){
          obj.MLE.old = obj.MLE
          
          z = (x.f - x.a.bar) / sqrt(n - 1)
          HZ = H %*% z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = lambda * z %*% t(z) %*% t(H) %*% HPHR.inv(lambda, HZ, R.inv, n)
          
          if(obj.MLE.old - obj.MLE < 1){
            break
          }else{
            x.a = x.f + K %*% d.f
            x.a.bar = apply(x.a, 1, mean)
          }
        }
        
        z = (x.f - x.a.bar) / sqrt(n - 1)
        # bw = cal_bandwidth_L96(z, GC)
        # P.f = dgTMat_L96(x.f, bw, GC)
        bw = Ccal_banding(as.matrix(t(z - apply(z, 1, mean))), p / 2)
        P.f = bandP(bw, p) * z %*% t(z)
        
        P.f.eigen = eigs(as.matrix(P.f), n)
        P.f.eigen.values = P.f.eigen$values
        P.f.eigen.vectors = P.f.eigen$vectors
        P.f.eigen.values[which(P.f.eigen.values < 1e-5)] = 0
        P.f = P.f.eigen.vectors %*% diag(P.f.eigen.values) %*% t(P.f.eigen.vectors)
        
        P.f.root = P.f.eigen.vectors %*% diag(P.f.eigen.values ^ 0.5)
        HZ = H %*% P.f.root
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = P.f %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        banding.bandwidth = bw
      }
      
      ### GC ####
      if(method == 'GC'){
        z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_GTE_bandwidth(t(z), dist_dim1_circ, GC, tol = 1)
        P.f = cal_taper_cov_L96(t(z), bw, dist_dim1_circ, GC)
        P.f.eigen = RSpectra::eigs(as.matrix(P.f), n)
        P.f.eigen.values = P.f.eigen$values
        P.f.eigen.vectors = P.f.eigen$vectors
        P.f.eigen.values[which(P.f.eigen.values < 1e-5)] = 0
        P.f = P.f.eigen.vectors %*% diag(P.f.eigen.values) %*% t(P.f.eigen.vectors)
        P.f.root = P.f.eigen.vectors %*% diag(P.f.eigen.values ^ 0.5)
        HZ = H %*% P.f.root
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = P.f %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = NA
        banding.bandwidth = bw
      }
      
      ### inflation then GC ####
      if(method == 'inflation then GC'){
        z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        bw = cal_GTE_bandwidth(t(z), dist_dim1_circ, GC, tol = 1)
        P.f = cal_taper_cov_L96(t(z), bw, dist_dim1_circ, GC)
        P.f.eigen = RSpectra::eigs(as.matrix(P.f), n)
        P.f.eigen.values = P.f.eigen$values
        P.f.eigen.vectors = P.f.eigen$vectors
        P.f.eigen.values[which(P.f.eigen.values < 1e-5)] = 0
        P.f = P.f.eigen.vectors %*% diag(P.f.eigen.values) %*% t(P.f.eigen.vectors)
        P.f.root = P.f.eigen.vectors %*% diag(P.f.eigen.values ^ 0.5)
        HZ = H %*% P.f.root
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = P.f %*% t(H) %*% HPHR.inv(1, HZ, R.inv, n)
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
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