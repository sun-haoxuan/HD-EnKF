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