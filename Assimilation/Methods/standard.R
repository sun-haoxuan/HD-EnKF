DAfun = function(forecast, DAinfo, i){
  x.f = forecast$x.f
  d.f = forecast$d.f
  H = forecast$H
  R = forecast$R
  R.inv = forecast$R.inv
  R.inv.root = forecast$R.inv.root
  x.f.bar = apply(x.f, 1, mean)
  d.f.bar = apply(d.f, 1, mean)

  Z = (x.f - x.f.bar) / sqrt(n - 1)
  HZ = H %*% Z
  D = svd(R.inv.root %*% HZ)$d
  obj.MLE = log(HPHR.det2(1, D, R, sd(D))) + length(D) * log(sd(D))+
    t(d.f.bar) %*% HPHR.inv(1, HZ, R.inv, ncol(Z)) %*% d.f.bar
  K = Z %*% t(HZ) %*% HPHR.inv(1, HZ, R.inv, ncol(Z))
  
  analysis = list(
    x.a = x.f + K %*% d.f,
    inflation.fator = 1,
    objective.likelihood = obj.MLE,
    iteration.number = NA,
    banding.bandwidth = NA
  )
  
  return(analysis)
}
