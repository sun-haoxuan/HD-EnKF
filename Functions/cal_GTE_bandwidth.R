risk_each_dist = function(x, dist.func, digits = 0, limit = Inf){
  p = ncol(x)
  n = nrow(x)
  
  tmp = list()
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      d_ij = round(dist.func(i, j, p), digits)
      if(d_ij <= limit){
        d_ij = as.character(d_ij)
        if(! d_ij %in% names(tmp)){
          tmp[[d_ij]] = c(0, 0)
        }
        # for(l in 1:n){
        #   for(m in setdiff(1:n, l)){
        #     tmp[[d_ij]][1] = tmp[[d_ij]][1] + (x[l, i] ^ 2) * (x[m, j] ^ 2) / n
        #     tmp[[d_ij]][2] = tmp[[d_ij]][2] + x[l, i] * x[l, j] * x[m, i] * x[m, j]
        #   }
        # }
        x_i_2 = x[, i] ^ 2
        x_j_2 = x[, j] ^ 2
        x_i_j = x[, i] * x[, j]
        xi2_xj2 = sum(x_i_2 * x_j_2)
        tmp[[d_ij]][1] = tmp[[d_ij]][1] + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
        tmp[[d_ij]][2] = tmp[[d_ij]][2] + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
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

cal_GTE_bandwidth = function(x, dist.func, taper.func, digits = 0, limit = Inf, tol = .Machine$double.eps^0.25){
  risk = risk_each_dist(x, dist.func, digits, limit)
  opt = optimise(function(k){risk_taper(risk, k, taper.func)}, c(0, max(as.numeric(names(risk)))), tol = tol)
  
  return(opt$minimum)
}