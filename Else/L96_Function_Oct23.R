RMSE = function(a, b){return(sqrt(mean((a - b) ^2 , na.rm = T)))}

plot.heatmap = function(x){corrplot::corrplot(as.matrix(x), method = 'color', is.corr = F)}

HPHR.det = function(lambda, D, R){
  return(det(R) * prod(lambda * D ^ 2 + 1))
}

HPHR.inv = function(lambda, HZ, R.inv, n){
  tmp = solve(diag(n) + lambda * t(HZ) %*% R.inv %*% HZ) %*% t(HZ) %*% R.inv
  return(R.inv - lambda * R.inv %*% HZ %*% tmp)
}

cal_bandwidth_L96 = function(x, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  # coor = matrix(c(cos(1:p * 2 * pi / p), sin(1:p * 2 * pi / p)), p, 2)
  coor = 1:p
  
  M_n = function(k){
    tmp = 0
    for(i in 1:(p - 1)){
      sigma_ii = sqrt(sum(x.scale[i, ] ^ 2))
      for(j in i:p){
        # grid.distance = sqrt(sum((coor[j, ] - coor[i, ]) ^ 2))
        grid.distance = min(coor[j] - coor[i], p - coor[j] + coor[i])
        if(grid.distance > k){
          tmp = tmp + sum(x.scale[i, ] * x.scale[j, ]) / p
        }else{
          tmp = tmp + sigma_ii * sqrt(sum(x.scale[j, ] ^ 2)) / p / n
        }
      }
    }
    return(tmp)
  }
  
  opt = optimise(M_n, c(0, p / 2), tol = 1)
  bw = opt$minimum
  return(bw)
}

cal_mapping_L96 = function(x, bw, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  
  # coor = matrix(c(cos(1:p * 2 * pi / p), sin(1:p * 2 * pi / p)), p, 2)
  coor = 1:p
  
  obj.svd = svd(x.scale)
  u = obj.svd$u
  d = obj.svd$d
  d[which(d < 1e-5)] = 1e-5
  x.scale.inv = t(x.scale) %*% u %*% diag(1 / d ^ 2) %*% t(u)
  
  # mapping.matrix = matrix(0, n, n)
  # for(k in 1:p){
  #   for(l in 1:p){
  #     grid.distance = sqrt(sum((coor[k, ] - coor[l, ]) ^ 2))
  #     if(grid.distance <= bw){
  #       Sigma2_kl = sum(x.scale[k, ] * x.scale[l, ])
  #       for(i in 1:n){
  #         for(j in 1:n){
  #           mapping.matrix[i, j] = mapping.matrix[i, j] + x.scale.inv[i, k] * x.scale.inv[j, l] * Sigma2_kl
  #         }
  #       }
  #     }
  #   }
  # }
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
      # grid.distance = sqrt(sum((coor[l, ] - coor[m, ]) ^ 2))
      grid.distance = min(coor[j] - coor[i], p - coor[j] + coor[i])
      if(grid.distance <= bw){
        Sigma2_lm = sum(x.scale[l, ] * x.scale[m, ])
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

RK4 = function(x, t, h, p, Fc){
  L96 = function(x, t){
    index1 = c(2:p, 1)
    index2 = c(p - 1, p, 1:(p - 2))
    index3 = c(p, 1:(p - 1))
    d = (x[index1] - x[index2]) * x[index3] - x[1:p] + Fc
    return(d)
  }
  k1 = h * L96(x, t)
  k2 = h * L96(x + 0.5 * k1, t + 0.5 * h)
  k3 = h * L96(x + 0.5 * k2, t + 0.5 * h)
  k4 = h * L96(x + k3, t + h)
  return(x + (k1 + 2 * k2 + 2 * k3 + k4) / 6)
}

L96_State = function(option){
  require(MASS)
  
  p = option$p
  S = option$S
  F.true = option$F.true
  F.set = option$F.set
  h = option$h
  sigma.tr = option$sigma.tr
  p.o = option$p.o
  S.o = option$S.o
  sigma = option$sigma
  rho = option$rho
  
  nsite = length(p.o)
  nobs = length(S.o)
  ti = seq(h, h * S, h)
  
  x.t = array(0, dim = c(p, length(ti)))
  x.t[, 1] = F.true; x.t[round(p / 2), 1] = F.true + 0.001
  for(i in 1:(S - 1)){
    x.t[, i + 1] = RK4(x.t[, i], ti[i], ti[i + 1] - ti[i], p, F.true)
  }
  
  H = matrix(0, nsite, p)
  for(k in 1:nsite){
    H[k, p.o[k]] = 1
  }
  
  R = matrix(0, nsite, nsite)
  for(i in 1:nsite){
    for(j in 1:nsite){
      R[i,j] = sigma ^ 2 * rho ^ min(abs(i - j), nsite - abs(i-j))
    }
  }
  
  y.o = array(NA, dim = c(p, nobs))
  y.o[p.o, 1:nobs] = x.t[p.o, S.o] +
    t(MASS::mvrnorm(nobs, rep(0, nsite), R))
  
  output = list(
    true = x.t,
    observation = y.o,
    param = list(p = p, S = S, F.true = F.true, h = h, 
                 p.o = p.o, S.o = S.o, sigma = sigma, rho = rho)
  )
  return(output)
}

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