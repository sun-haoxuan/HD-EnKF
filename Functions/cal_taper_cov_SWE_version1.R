cal_taper_cov_SWE_v1 = function(x, coor, taper.func, digits = 1, tol = 0.1){
  n = nrow(x)
  p = ncol(x)
  ncoor = nrow(coor)
  
  # t1 = Sys.time()
  mat.dist = diag(0, ncoor)
  for(i in 1:(ncoor - 1)){
    for(j in (i + 1):ncoor){
      d_ij = round(sqrt(sum((coor[i, ] - coor[j, ]) ^ 2)), digits)
      mat.dist[i, j] = mat.dist[j, i] = d_ij
    }
  }

  risk = list()
  dist.all = unique(as.vector(mat.dist))
  for(d_ij in dist.all){
    list.grid = which(mat.dist == d_ij, arr.ind = TRUE)
    tmp1 = tmp2 = 0
    for(grd in 1:nrow(list.grid)){
      grid.x = list.grid[grd, 1]
      grid.y = list.grid[grd, 2]
      for(nvar1 in 0:2){
        for(nvar2 in 0:2){
          i = ncoor * nvar1 + grid.x
          j = ncoor * nvar2 + grid.y
          # print(paste(i, j))
          x_i_2 = x[, i] ^ 2
          x_j_2 = x[, j] ^ 2
          x_i_j = x[, i] * x[, j]
          xi2_xj2 = sum(x_i_2 * x_j_2)
          tmp1 = tmp1 + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
          tmp2 = tmp2 + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
        }
      }
    }
    risk[[as.character(d_ij)]] = c(tmp1, tmp2)
  }
  # t2 = Sys.time()
  # print(t2 - t1)
  # save(risk, file = 'test_risk_newcode2.Rdata')
  
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
  
  taper.optimise = optimise(risk_taper, c(0, max(as.numeric(names(risk)))), tol = tol)
  taper.bandwidth = taper.optimise$minimum

  mat.dist.taper.0 = diag(1, ncoor)
  for(i in 1:(ncoor - 1)){
    for(j in (i + 1):ncoor){
      mat.dist.taper.0[i, j] = mat.dist.taper.0[j, i] = taper.func(mat.dist[i, j] / taper.bandwidth)
    }
  }
  mat.dist.taper = matrix(NA, 3 * ncoor, 3 * ncoor)
  for(i in 0:2){
    for(j in 0:2){
      mat.dist.taper[ncoor * i + 1:ncoor, ncoor * j + 1:ncoor] = mat.dist.taper.0
    }
  }
  
  Sigma.hat = mat.dist.taper * (t(x) %*% x)
  Sigma.hat.eigen = RSpectra::eigs(Sigma.hat, n)
  Sigma.hat.eigen.values = Sigma.hat.eigen$values
  Sigma.hat.eigen.vectors = Sigma.hat.eigen$vectors
  Sigma.hat.eigen.values[which(Sigma.hat.eigen.values < 1e-5)] = 0
  Sigma.tilde.root = Sigma.hat.eigen.vectors %*% diag(Sigma.hat.eigen.values ^ 0.5)
  
  output = list(
    bandwidth = taper.bandwidth,
    Z = Sigma.tilde.root
  )
  
  return(output)
}