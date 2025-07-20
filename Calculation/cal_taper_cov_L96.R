banding = BL = function(z){
  if(z <= 1){
    return(1)
  }else{
    return(0)
  }
}

tapering = CZZ = function(z){
  if(z <= 0.5){
    return(1)
  }else if(z > 0.5 & z <=1){
    return(2 - 2 * z)  
  }else{
    return(0)
  }
}

GC = function(z){
  z = z * 2
  if(z >= 0 & z <= 1){
    return(1 - 5 * z ^ 2 / 3 + 5 * z ^ 3 / 8 + z ^ 4 / 2 - z ^ 5 / 4)
  }else if(z > 1 & z <= 2){
    return(- 2 / 3 / z + 4 - 5 * z + 5 * z ^ 2 / 3 + 5 * z ^ 3 / 8 - z ^ 4 / 2 + z ^ 5 / 12)
  }else{
    return(0)
  }
}

cal_taper_bandwidth = function(x, mat.dist, taper.func, interval, tol = .Machine$double.eps^0.25){
  n = nrow(x)
  p = ncol(x)
  
  risk = list()
  dist.all = unique(as.vector(mat.dist))
  for(d_ij in dist.all){
    list.grid = which(mat.dist == d_ij, arr.ind = TRUE)
    tmp1 = tmp2 = 0
    for(grd in 1:nrow(list.grid)){
      i = list.grid[grd, 1]
      j = list.grid[grd, 2]
      x_i_2 = x[, i] ^ 2
      x_j_2 = x[, j] ^ 2
      x_i_j = x[, i] * x[, j]
      xi2_xj2 = sum(x_i_2 * x_j_2)
      tmp1 = tmp1 + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
      tmp2 = tmp2 + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
    }
    risk[[as.character(d_ij)]] = c(tmp1, tmp2)
  }
  
  risk_taper = function(k){
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
  
  taper.optimise = optimise(risk_taper, c(interval[1], interval[2]), tol = tol)
  taper.bandwidth = taper.optimise$minimum

  return(taper.bandwidth)
}

cal_taper_Z = function(x, taper.bandwidth, mat.dist, taper.func, ndim = min(n, p)){
  n = nrow(x)
  p = ncol(x)
  
  Sigma.hat = t(x) %*% x
  mat.dist.taper = diag(1, p)
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      mat.dist.taper[i, j] = mat.dist.taper[j, i] = taper.func(mat.dist[i, j] / taper.bandwidth)
    }
  }
  Sigma.hat = Sigma.hat * mat.dist.taper
  Sigma.hat.eigen = RSpectra::eigs(Sigma.hat, ndim)
  Sigma.hat.eigen.values = Sigma.hat.eigen$values
  Sigma.hat.eigen.vectors = Sigma.hat.eigen$vectors
  Sigma.hat.eigen.values[which(Sigma.hat.eigen.values < 1e-5)] = 1e-5
  Sigma.tilde.root = Sigma.hat.eigen.vectors %*% diag(Sigma.hat.eigen.values ^ 0.5)
  
  return(Sigma.tilde.root)
}