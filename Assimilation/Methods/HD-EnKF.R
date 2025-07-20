DAfun = function(forecast, DAinfo, i){
  x.f = forecast$x.f
  d.f = forecast$d.f
  H = forecast$H
  R = forecast$R
  R.inv = forecast$R.inv
  R.inv.root = forecast$R.inv.root
  x.f.bar = apply(x.f, 1, mean)
  d.f.bar = apply(d.f, 1, mean)
  
  taper = DAinfo$taper
  interval = DAinfo$interval
  num_eigen = DAinfo$num_eigen
  distMat = DAinfo$distMat
  
  Z = (x.f - x.f.bar) / sqrt(n - 1)
  bw = cal_taper_bandwidth(t(Z), distMat, taper, interval)
  Z = cal_taper_Z(t(Z), bw, distMat, taper, ndim = num_eigen)
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
  
  for(k in 1:10){
    x.a = x.f + K %*% d.f
    x.a.bar = apply(x.a, 1, mean)
    obj.MLE.old = obj.MLE

    Z = (x.f - x.a.bar) / sqrt(n - 1)
    Z = cal_taper_Z(t(Z), bw, distMat, taper, ndim = num_eigen)
    HZ = H %*% Z
    D = svd(R.inv.root %*% HZ)$d
    log.likelihood = function(lambda){
      return(log(HPHR.det2(lambda, D, R, sd(D))) + length(D) * log(sd(D))+
               t(d.f.bar) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z)) %*% d.f.bar)
    }
    likelihood.optimize = optimize(log.likelihood, c(0.1, 10))
    lambda = likelihood.optimize$minimum
    obj.MLE = likelihood.optimize$objective
    
    if(obj.MLE.old - obj.MLE < 1){
      break
    }else{
      K = lambda * Z %*% t(HZ) %*% HPHR.inv(lambda, HZ, R.inv, ncol(Z))
    }
  }

  analysis = list(
    x.a = x.f + K %*% d.f,
    inflation.fator = lambda0,
    objective.likelihood = obj.MLE.old,
    iteration.number = k,
    banding.bandwidth = bw
  )
  
  return(analysis)
}