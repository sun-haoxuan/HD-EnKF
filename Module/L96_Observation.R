L96_Observation = function(x.t, R, option){
  require(MASS)
  
  set.seed(option$boot)
  
  p = option$p
  S = option$S
  S.seq = option$S.seq
  
  p.o = sort(sample(p, q))
  S.o = seq(S.seq, S, S.seq)
  
  set.seed(b)
  nsite = length(p.o)
  nobs = length(S.o)

  y.o = array(NA, dim = c(p, nobs))
  y.o[p.o, 1:nobs] = x.t[p.o, S.o] +
    t(MASS::mvrnorm(nobs, rep(0, nsite), R))
  
  return(y.o)
}
