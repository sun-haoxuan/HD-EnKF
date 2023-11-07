cal_taper_cov_L96 = function(x, bw, dist.func, taper.func, scale = FALSE){
  n = nrow(x)
  p = ncol(x)
  
  if(scale){
    x = (x - apply(x, 2, mean)) / sqrt(n - 1)
  }
  
  taper.cov = matrix(0, p, p)
  for(i in 1:p){
    taper.cov[i, i] = sum(x[, i] ^ 2)
  }
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      grid.distance = dist.func(i, j, p)
      if(grid.distance <= bw){
        taper.cov[i, j] = taper.cov[j,  i] = sum(x[, i] * x[, j]) *  taper.func(grid.distance / bw)
      }
    }
  }
  
  return(taper.cov)
}