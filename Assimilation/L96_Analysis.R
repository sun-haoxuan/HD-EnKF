L96_Analysis = function(option){
  
  require(MASS)
  
  set.seed(option$boot)
  
  ## Initialize ####
  p = option$p
  q = option$q
  S = option$S
  S.seq = option$S.seq
  h = option$h
  
  ti = seq(h, h * S, h)
  p.o = sort(sample(p, q))
  S.o = seq(S.seq, S, S.seq)
  nobs = length(S.o)
  
  x.t = option$state$x.t
  y.o = option$state$y.o
  H = option$state$H
  R = option$state$R
  R.inv = option$state$R.inv
  R.inv.root = option$state$R.inv.root
  
  F.set = option$DA$F.set
  n = option$DA$n
  sigma.in = option$DA$sigma.in
  DAfun = option$DA$DAfun
  DAinfo <<- option$DA$DAinfo
  
  ## Output list ####
  cv = 0
  num = 0
  k = 0
  assimilate.errors = assimilate.spread = rep(NA, S)
  x.en.bar = matrix(NA, p, S)
  x.en.temp = x.t[, 1] + t(MASS::mvrnorm(n, rep(0, p), diag(p) * sigma.in))
  x.en.bar[, 1] = apply(x.en.temp, 1, mean)
  
  ## Analyse State ####
  assimilate.errors[1] = RMSE(x.t[, 1], x.en.bar[, 1])
  assimilate.spread[1] = cal_spread(x.en.temp)
  vec.inflation.factor = rep(NA, nobs)
  vec.objective.likelihood = rep(NA, nobs) 
  vec.iteration.number = rep(NA, nobs)
  vec.banding.bandwidth = rep(NA, nobs)
  
  ## Assimilation Methods ####
  forecast = list(
    x.f = NULL,
    d.f = NULL,
    H = H,
    R = R,
    R.inv = R.inv,
    R.inv.root = R.inv.root
  )
  
  for(i in 1:(S - 1)){
    ## Update by RK4 
    for(j in 1:n){
      x.en.temp[, j] = RK4(x.en.temp[, j], ti[i], ti[i + 1] - ti[i], p, F.set)
    }
    x.en.bar[, i + 1] = apply(x.en.temp, 1, mean)
    
    # Analyse if there exist observations
    if((i + 1) %in% S.o){
      num = num + 1
      
      forecast$x.f = x.en.temp
      forecast$d.f = y.o[p.o, num]  - H %*% x.en.temp + t(MASS::mvrnorm(n, rep(0, q), R))
      analysis = DAfun(forecast, DAinfo, i)
      
      x.en.temp = analysis$x.a
      x.en.bar[, i + 1] = apply(x.en.temp, 1, mean)
      vec.inflation.factor[num] = analysis$inflation.fator
      vec.objective.likelihood[num] = analysis$objective.likelihood
      vec.iteration.number[num] = analysis$iteration.number
      vec.banding.bandwidth[num] = analysis$banding.bandwidth
    }
    
    assimilate.errors[i + 1] = RMSE(x.t[, i + 1], x.en.bar[, i + 1])
    assimilate.spread[i + 1] = cal_spread(x.en.temp)
    print(paste(i + 1, round(assimilate.errors[i + 1], 2), round(assimilate.spread[i + 1], 2)))
    
    # Large error lead to divergence
    if(assimilate.errors[i + 1] > 10 | is.nan(assimilate.errors[i + 1])){
      cv = 1
      break
    }
  }
  print(paste('RMSE', round(sqrt(mean(assimilate.errors[1001:2000] ^ 2)), 3)))
  
  output = list(
    analyse = x.en.bar, 
    converge = cv, 
    errors = assimilate.errors,
    spread = assimilate.spread,
    inflation = vec.inflation.factor,
    likelihood = vec.objective.likelihood,
    iteration = vec.iteration.number,
    bandwidth = vec.banding.bandwidth
  )
  
  return(output)
}