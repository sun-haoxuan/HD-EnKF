risk_each_dist = function(x, dist.func){
  n = nrow(x)
  p = ncol(x)
  
  tmp = list()
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      d_ij = as.character(dist.func(i, j, p))
      if(! d_ij %in% names(tmp)){
        tmp[[d_ij]] = c(0, 0)
      }
      for(l in 1:n){
        for(m in setdiff(1:n, l)){
          tmp[[d_ij]][1] = tmp[[d_ij]][1] + (x[l, i] ^ 2) * (x[m, j] ^ 2) / n
          tmp[[d_ij]][2] = tmp[[d_ij]][2] + x[l, i] * x[l, j] * x[m, i] * x[m, j]
        }
      }
    }
  }
  return(tmp)
}

risk_taper = function(risk, k, taper.func){
  list.dist = names(risk)[order(as.numeric(names(risk)))]
  tmp = 0
  id_dist = max(which(as.numeric(list.dist) <= k))
  if(!is.infinite(id_dist)){
    for(i in 1:id_dist){
      g_ij = taper.func(as.numeric(list.dist)[i] / k)
      tmp = tmp + g_ij ^ 2 * risk[[list.dist[i]]][1]  + (g_ij ^ 2 - 2 * g_ij) * risk[[list.dist[i]]][2]
    }
  }
  
  return(tmp)
}

cal_GTE_bandwidth = function(x, dist.func, taper.func, tol = .Machine$double.eps^0.25){
  risk = risk_each_dist(x, dist.func)
  opt = optimise(function(k){risk_taper(risk, k, taper.func)}, c(0, max(as.numeric(names(risk)))), tol = tol)
  
  return(opt$minimum)
}