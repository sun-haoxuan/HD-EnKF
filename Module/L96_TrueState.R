L96_TrueState = function(option){
  p = option$p
  S = option$S
  F.true = option$F.true
  h = option$h

  ti = seq(h, h * S, h)
  
  ## true state
  x.t = array(0, dim = c(p, length(ti)))
  x.t[, 1] = F.true; x.t[round(p / 2), 1] = F.true + 0.001
  for(i in 1:(S - 1)){
    x.t[, i + 1] = RK4(x.t[, i], ti[i], ti[i + 1] - ti[i], p, F.true)
  }
  
  return(x.t)
}