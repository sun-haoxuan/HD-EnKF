L96_Analysis = function(state, option){
  
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
  
  taper = option$taper
  is.oracle = option$is.oracle
  is.inflation = option$is.inflation
  is.iteration = option$is.iteration
  interval = option$interval
  neigen = option$neigen
  
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
  vec.number.eigen = rep(NA, nobs)
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
      
      ### Assimilate ####
      lambda0 = obj.MLE = bw = k = NA
      Z = (x.f - x.f.bar) / sqrt(n - 1)
      bw = cal_taper_bandwidth(t(Z), mat.dist, taper, interval)
      Z = cal_taper_Z(t(Z), bw, mat.dist, taper, ndim = p)
      
      A1 = svd(Z)$d ^ 2
      A2 = cumsum(A1) / sum(A1)
      A3 = which(A2 > 0.9)[1]
      Z = cal_taper_Z(t(Z), bw, mat.dist, taper, ndim = A3)
      
      HZ = H %*% Z
      D = svd(R.inv.root %*% HZ)$d
      log.likelihood = function(lambda){
        return(log(HPHR.det2(lambda, D, R, sd(D) * lambda)) + length(D) * log(sd(D) * lambda)+
                 t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
      }
      likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
      lambda0 = lambda = likelihood.optimize$minimum
      obj.MLE = likelihood.optimize$objective
      K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
      
      x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
      number.eigen = A3
      
      k = 0
      for(k in 1:10){
        x.a = x.f + K %*% d.f
        x.a.bar = apply(x.a, 1, mean)
        obj.MLE.old = obj.MLE
        K.old = K
        
        Z = (x.f - x.a.bar) / sqrt(n - 1)
        Z = cal_taper_Z(t(Z), bw, mat.dist, taper, ndim = p)
        
        A1 = svd(Z)$d ^ 2
        A2 = cumsum(A1) / sum(A1)
        A3 = which(A2 > 0.9)[1]
        Z = cal_taper_Z(t(Z), bw, mat.dist, taper, ndim = A3)
        
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det2(lambda, D, R, sd(D) * lambda)) + length(D) * log(sd(D) * lambda)+
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
      
      ### Update ####
      x.en.temp = x.f + K %*% d.f
      x.en.bar[, i + 1] = apply(x.en.temp, 1, mean)
      vec.inflation.factor[num] = lambda0
      vec.objective.likelihood[num] = obj.MLE
      vec.iteration.number[num] = k
      vec.banding.bandwidth[num] = bw
      vec.number.eigen[num] = number.eigen
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