#   nx = 50
#   ny = 31
#   coor = data.frame(
#     x = rep(1:nx, times = ny),
#     y = rep(1:ny, each = nx)
#   )
#   ncoor = nrow(coor)
#   mat.dist = diag(0, ncoor)
#   for(i in 1:(ncoor - 1)){
#     for(j in (i + 1):ncoor){
#       d_ij = round(sqrt(sum((coor[i, ] - coor[j, ]) ^ 2)), 1)
#       mat.dist[i, j] = mat.dist[j, i] = d_ij
#     }
#   }
#   save(nx, ny, coor, mat.dist, file = 'Output/matrix_distance.Rdata')
load('Output/matrix_distance.Rdata')

banding = function(z){
  if(z <= 1){
    return(1)
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

cal_taper_bandwidth_block = function(x, mat.dist, taper.func = GC, interval, tol = .Machine$double.eps^0.25){
  n = nrow(x)
  p = ncol(x)
  
  ncoor = nrow(mat.dist)
  
  vec.taper.bandwidth = c()
  for(nvar1 in 0:2){
    for(nvar2 in nvar1:2){
      print(paste(nvar1, nvar2))
      
      risk = list()
      dist.all = unique(as.vector(mat.dist))
      for(d_ij in dist.all){
        list.grid = which(mat.dist == d_ij, arr.ind = TRUE)
        tmp1 = tmp2 = 0
        for(grd in 1:nrow(list.grid)){
          i = ncoor * nvar1 + list.grid[grd, 1]
          j = ncoor * nvar2 + list.grid[grd, 2]
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
      vec.taper.bandwidth = c(vec.taper.bandwidth, taper.optimise$minimum)
    }
  }
  print(vec.taper.bandwidth)
  
  return(vec.taper.bandwidth)
}

cal_taper_bandwidth_total = function(x, mat.dist, taper.func = GC, interval, tol = .Machine$double.eps^0.25){
  n = nrow(x)
  p = ncol(x)
  
  ncoor = nrow(mat.dist)
  
  risk = list()
  dist.all = unique(as.vector(mat.dist))
  for(d_ij in dist.all){
    # print(d_ij)
    list.grid = which(mat.dist == d_ij, arr.ind = TRUE)
    tmp1 = tmp2 = 0
    for(grd in 1:nrow(list.grid)){
      for(nvar1 in 0:2){
        for(nvar2 in nvar1:2){
          i = ncoor * nvar1 + list.grid[grd, 1]
          j = ncoor * nvar2 + list.grid[grd, 2]
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
  print(taper.bandwidth)
  
  return(taper.bandwidth)
}

cal_taper_Z_block = function(x, vec.taper.bandwidth, taper.func = GC, tol = .Machine$double.eps^0.25){
  n = nrow(x)
  p = ncol(x)
  ncoor = nrow(mat.dist)
  
  Sigma.hat = t(x) %*% x
  n.bw = 0
  for(nvar1 in 0:2){
    for(nvar2 in nvar1:2){
      # print(paste(nvar1, nvar2))
      n.bw = n.bw + 1
      mat.dist.taper = diag(1, ncoor)
      for(i in 1:(ncoor - 1)){
        for(j in (i + 1):ncoor){
          mat.dist.taper[i, j] = mat.dist.taper[j, i] = taper.func(mat.dist[i, j] / vec.taper.bandwidth[n.bw])
        }
      }
      Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor]
      Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor]
    }
  }
  Sigma.hat.eigen = RSpectra::eigs(Sigma.hat, n)
  Sigma.hat.eigen.values = Sigma.hat.eigen$values
  Sigma.hat.eigen.vectors = Sigma.hat.eigen$vectors
  Sigma.hat.eigen.values[which(Sigma.hat.eigen.values < 1e-5)] = 0
  Sigma.tilde.root = Sigma.hat.eigen.vectors %*% diag(Sigma.hat.eigen.values ^ 0.5)
  
  return(Sigma.tilde.root)
}

cal_taper_Z_total = function(x, taper.bandwidth, taper.func = GC, tol = .Machine$double.eps^0.25){
  n = nrow(x)
  p = ncol(x)
  ncoor = nrow(mat.dist)
  
  Sigma.hat = t(x) %*% x
  mat.dist.taper = diag(1, ncoor)
  for(i in 1:(ncoor - 1)){
    for(j in (i + 1):ncoor){
      mat.dist.taper[i, j] = mat.dist.taper[j, i] = taper.func(mat.dist[i, j] / taper.bandwidth)
    }
  }
  for(nvar1 in 0:2){
    Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar1 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar1 + 1:ncoor]
  }
  for(nvar1 in 0:2){
    for(nvar2 in nvar1:2){
      Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor]
      Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor]
    }
  }
  Sigma.hat.eigen = RSpectra::eigs(Sigma.hat, n)
  Sigma.hat.eigen.values = Sigma.hat.eigen$values
  Sigma.hat.eigen.vectors = Sigma.hat.eigen$vectors
  Sigma.hat.eigen.values[which(Sigma.hat.eigen.values < 1e-5)] = 0
  Sigma.tilde.root = Sigma.hat.eigen.vectors %*% diag(Sigma.hat.eigen.values ^ 0.5)
  
  return(Sigma.tilde.root)
}

cal_taper_cov_SWE_block = function(x, mat.dist, taper.func = GC, tol = .Machine$double.eps^0.25){
  n = nrow(x)
  p = ncol(x)
  ncoor = nrow(mat.dist)
  
  vec.taper.bandwidth = c()
  for(nvar1 in 0:2){
    for(nvar2 in nvar1:2){
      # print(paste(nvar1, nvar2))
      
      risk = list()
      dist.all = unique(as.vector(mat.dist))
      for(d_ij in dist.all){
        list.grid = which(mat.dist == d_ij, arr.ind = TRUE)
        tmp1 = tmp2 = 0
        for(grd in 1:nrow(list.grid)){
          i = ncoor * nvar1 + list.grid[grd, 1]
          j = ncoor * nvar2 + list.grid[grd, 2]
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
      
      taper.optimise = optimise(risk_taper, c(0, max(as.numeric(names(risk)))), tol = tol)
      vec.taper.bandwidth = c(vec.taper.bandwidth, taper.optimise$minimum)
    }
  }
  
  Sigma.hat = t(x) %*% x
  n.bw = 0
  for(nvar1 in 0:2){
    for(nvar2 in nvar1:2){
      print(paste(nvar1, nvar2))
      n.bw = n.bw + 1
      mat.dist.taper = diag(1, ncoor)
      for(i in 1:(ncoor - 1)){
        for(j in (i + 1):ncoor){
          mat.dist.taper[i, j] = mat.dist.taper[j, i] = taper.func(mat.dist[i, j] / vec.taper.bandwidth[n.bw])
        }
      }
      Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor]
      Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor]
    }
  }
  Sigma.hat.eigen = RSpectra::eigs(Sigma.hat, n)
  Sigma.hat.eigen.values = Sigma.hat.eigen$values
  Sigma.hat.eigen.vectors = Sigma.hat.eigen$vectors
  Sigma.hat.eigen.values[which(Sigma.hat.eigen.values < 1e-5)] = 0
  Sigma.tilde.root = Sigma.hat.eigen.vectors %*% diag(Sigma.hat.eigen.values ^ 0.5)
  
  output = list(
    bandwidth = vec.taper.bandwidth,
    Z = Sigma.tilde.root
  )
  
  return(output)
}

cal_taper_cov_SWE_total = function(x, mat.dist, taper.func = GC, tol = .Machine$double.eps^0.25){
  n = nrow(x)
  p = ncol(x)
  
  ncoor = nrow(mat.dist)
  
  risk = list()
  dist.all = unique(as.vector(mat.dist))
  for(d_ij in dist.all){
    print(d_ij)
    list.grid = which(mat.dist == d_ij, arr.ind = TRUE)
    tmp1 = tmp2 = 0
    for(grd in 1:nrow(list.grid)){
      for(nvar1 in 0:2){
        for(nvar2 in nvar1:2){
          i = ncoor * nvar1 + list.grid[grd, 1]
          j = ncoor * nvar2 + list.grid[grd, 2]
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
  
  
  Sigma.hat = t(x) %*% x
  mat.dist.taper = diag(1, ncoor)
  for(i in 1:(ncoor - 1)){
    for(j in (i + 1):ncoor){
      mat.dist.taper[i, j] = mat.dist.taper[j, i] = taper.func(mat.dist[i, j] / taper.bandwidth)
    }
  }
  for(nvar1 in 0:2){
    Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar1 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar1 + 1:ncoor]
  }
  for(nvar1 in 0:2){
    for(nvar2 in nvar1:2){
      Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar1 + 1:ncoor, ncoor * nvar2 + 1:ncoor]
      Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor] = mat.dist.taper * Sigma.hat[ncoor * nvar2 + 1:ncoor, ncoor * nvar1 + 1:ncoor]
    }
  }
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

# t1 = Sys.time(); A = cal_taper_cov_SWE_total(x, mat.dist); t2 = Sys.time(); print(t2 - t1)
