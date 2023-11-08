SWE_Analyse = function(state, opt){
  
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
  boot = opt$boot
  
  ## Other parameters
  ngrid <- nx * ny
  dx <- L / nx
  dy <- D / (ny - 1)
  nobs <- (S.o.end - S.o.start) %/% S.o.seq + 1 # forecast 1h, then assimilate
  ndim <- ngrid * 3
  nid <- length(id)
  xobs <- ny * nid * 3 # observe h only
  
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
  vec.taper.bandwidth = matrix(NA, 6, nobs)
  assimilate.errors = matrix(NA, 4, S)
  assimilate.errors[, 1] = c(RMSE(x.t[, 1], x.en.bar[, 1]),
                           RMSE(x.t[1:ncoor, 1], x.en.bar[1:ncoor, 1]), 
                           RMSE(x.t[1:ncoor + ncoor, 1], x.en.bar[1:ncoor + ncoor, 1]), 
                           RMSE(x.t[1:ncoor + ncoor * 2, 1], x.en.bar[1:ncoor + ncoor * 2, 1]))
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
      
      ### standard ####
      if(method == 'standard'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        mulo = prod(1 + svd(R.inv.root %*% HZ)$d ^ 2)
        obj.MLE = logdetR + log(mulo) + t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = NA
      }
      
      ### inflation & iteration ####
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
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = NA
      }
      
      ### taper-inflation ####
      if(method == 'taper-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, GC)
        Z = cal_taper_Z_total(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### taper-iteration ####
      if(method == 'taper-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, GC)
        Z = cal_taper_Z_total(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_total(t(Z), bw, GC)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          obj.MLE = log(HPHR.det(1, D, R)) +
            t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### taper-iteration(withinflation) ####
      if(method == 'taper-iteration(withinflation)'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, GC)
        Z = cal_taper_Z_total(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_total(t(Z), bw, GC)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### taper-inflation-iteration ####
      if(method == 'taper-inflation-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, GC)
        Z = cal_taper_Z_total(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_total(t(Z), bw, GC)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda0 = lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### inflation-taper ####
      if(method == 'inflation-taper'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum

        bw = cal_taper_bandwidth_total(t(Z), mat.dist, GC)
        Z = cal_taper_Z_total(sqrt(lambda) * t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar

        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### inflation-taper-iteration ####
      if(method == 'inflation-taper-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, GC)
        Z = cal_taper_Z_total(sqrt(lambda0) * t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_total(sqrt(lambda0) * t(Z), bw, GC)
          D = svd(R.inv.root %*% HZ)$d
          obj.MLE = log(HPHR.det(1, D, R)) +
            t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### taper(B)-inflation ####
      if(method == 'taper(B)-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, GC)
        Z = cal_taper_Z_block(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### taper(B)-iteration ####
      if(method == 'taper(B)-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, GC)
        Z = cal_taper_Z_block(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_block(t(Z), bw, GC)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          obj.MLE = log(HPHR.det(1, D, R)) +
            t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### taper(B)-iteration(withinflation) ####
      if(method == 'taper(B)-iteration(withinflation)'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, GC)
        Z = cal_taper_Z_block(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_block(t(Z), bw, GC)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### taper(B)-inflation-iteration ####
      if(method == 'taper(B)-inflation-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, GC)
        Z = cal_taper_Z_block(t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_block(t(Z), bw, GC)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda0 = lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        x.f = sqrt(lambda0) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda0
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### inflation-taper(B) ####
      if(method == 'inflation-taper(B)'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, GC)
        Z = cal_taper_Z_block(sqrt(lambda) * t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### inflation-taper(B)-iteration ####
      if(method == 'inflation-taper(B)-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda0 = lambda = likelihood.optimize$minimum
        
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, GC)
        Z = cal_taper_Z_block(sqrt(lambda0) * t(Z), bw, GC)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_block(sqrt(lambda0) * t(Z), bw, GC)
          D = svd(R.inv.root %*% HZ)$d
          obj.MLE = log(HPHR.det(1, D, R)) +
            t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### band-inflation ####
      if(method == 'band-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_total(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### band-iteration ####
      if(method == 'band-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_total(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_total(t(Z), bw, banding)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          obj.MLE = log(HPHR.det(1, D, R)) +
            t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### band-iterwithinfl ####
      if(method == 'band-iterwithinfl'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_total(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_total(t(Z), bw, banding)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### inflation-band ####
      if(method == 'inflation-band'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_total(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### inflation-bandwithinfl ####
      if(method == 'inflation-bandwithinfl'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        
        bw = cal_taper_bandwidth_total(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_total(sqrt(lambda) * t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### bandblock-inflation ####
      if(method == 'bandblock-inflation'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_block(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        obj.MLE = likelihood.optimize$objective
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### bandblock-iteration ####
      if(method == 'bandblock-iteration'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_block(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_block(t(Z), bw, banding)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          obj.MLE = log(HPHR.det(1, D, R)) +
            t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
          K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### bandblock-iterwithinfl ####
      if(method == 'bandblock-iterwithinfl'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_block(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        for(k in 1:10){
          x.a = x.f + K %*% d.f
          x.a.bar = apply(x.a, 1, mean)
          obj.MLE.old = obj.MLE
          K.old = K
          
          Z = (x.f - x.a.bar) / sqrt(n - 1)
          Z = cal_taper_Z_block(t(Z), bw, banding)
          HZ = H %*% Z
          D = svd(R.inv.root %*% HZ)$d
          log.likelihood = function(lambda){
            return(log(HPHR.det(lambda, D, R)) +
                     t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
          }
          likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
          lambda = likelihood.optimize$minimum
          obj.MLE = likelihood.optimize$objective
          K = Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, n)
          
          print(paste(k, obj.MLE.old, obj.MLE))
          if(obj.MLE.old - obj.MLE < 1){
            break
          }
        }
        obj.MLE = obj.MLE
        K = K.old
        
        inflation.fator = NA
        objective.likelihood = obj.MLE
        iteration.number = k
        taper.bandwidth = bw
      }
      
      ### inflation-bandblock ####
      if(method == 'inflation-bandblock'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_block(t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }
      
      ### inflation-bandblockwithinfl ####
      if(method == 'inflation-bandblockwithinfl'){
        Z = (x.f - x.f.bar) / sqrt(n - 1)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        log.likelihood = function(lambda){
          return(log(HPHR.det(lambda, D, R)) +
                   t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, n) %*% d.f.bar)
        }
        likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
        lambda = likelihood.optimize$minimum
        
        bw = cal_taper_bandwidth_block(t(Z), mat.dist, banding, c(0, max(mat.dist)))
        Z = cal_taper_Z_block(sqrt(lambda) * t(Z), bw, banding)
        HZ = H %*% Z
        D = svd(R.inv.root %*% HZ)$d
        obj.MLE = log(HPHR.det(1, D, R)) +
          t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, n) %*% d.f.bar
        K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, n)
        
        x.f = sqrt(lambda) * (x.f - x.f.bar) + x.f.bar
        
        inflation.fator = lambda
        objective.likelihood = obj.MLE
        iteration.number = NA
        taper.bandwidth = bw
      }

      ### Update ####
      x.en.temp = x.f + K %*% d.f
      x.en.bar[, t + 1] = apply(x.en.temp, 1, mean)
      vec.inflation.factor[num] = inflation.fator
      vec.objective.likelihood[num] = objective.likelihood
      vec.iteration.number[num] = iteration.number
      vec.taper.bandwidth[, num] = taper.bandwidth
      
      print(paste('Boot', boot, 'Step', t, 'error', round(RMSE(x.t[, t + 1], x.en.bar[, t + 1]), 3)))
    }
    assimilate.errors[, t + 1] = c(RMSE(x.t[, t + 1], x.en.bar[, t + 1]),
                               RMSE(x.t[1:ncoor, t + 1], x.en.bar[1:ncoor, t + 1]), 
                               RMSE(x.t[1:ncoor + ncoor, t + 1], x.en.bar[1:ncoor + ncoor, t + 1]), 
                               RMSE(x.t[1:ncoor + ncoor * 2, t + 1], x.en.bar[1:ncoor + ncoor * 2, t + 1]))
    
    if(assimilate.errors[1, t + 1] > 20 | is.na(assimilate.errors[1, t + 1])){
      cv = 1
      break
    }
  }
  
  output = list(
    analyse = x.en.bar,
    converge = cv, 
    error = assimilate.errors, 
    objective = vec.objective.likelihood,
    iteration = vec.iteration.number,
    inflation = vec.inflation.factor,
    bandwidth = vec.taper.bandwidth,
    option = opt 
    )
  
  return(output)
}