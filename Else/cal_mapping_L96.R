cal_mapping_L96 = function(x, bw, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  
  # coor = matrix(c(cos(1:p * 2 * pi / p), sin(1:p * 2 * pi / p)), p, 2)
  coor = 1:p
  
  obj.svd = svd(x.scale)
  u = obj.svd$u
  d = obj.svd$d
  d[which(d < 1e-5)] = 1e-5
  x.scale.inv = t(x.scale) %*% u %*% diag(1 / d ^ 2) %*% t(u)
  
  # mapping.matrix = matrix(0, n, n)
  # for(k in 1:p){
  #   for(l in 1:p){
  #     grid.distance = sqrt(sum((coor[k, ] - coor[l, ]) ^ 2))
  #     if(grid.distance <= bw){
  #       Sigma2_kl = sum(x.scale[k, ] * x.scale[l, ])
  #       for(i in 1:n){
  #         for(j in 1:n){
  #           mapping.matrix[i, j] = mapping.matrix[i, j] + x.scale.inv[i, k] * x.scale.inv[j, l] * Sigma2_kl
  #         }
  #       }
  #     }
  #   }
  # }
  mapping.matrix = matrix(0, n, n)
  for(l in 1:p){
    Sigma2_ll = sum(x.scale[l, ] ^ 2)
    for(i in 1:n){
      mapping.matrix[i, i] = mapping.matrix[i, i] + x.scale.inv[i, l] ^ 2 * Sigma2_ll
    }
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
        mapping.matrix[i, j] = mapping.matrix[j, i] = mapping.matrix[i, j] + x.scale.inv[i, l] * x.scale.inv[j, l] * Sigma2_ll
      }
    }
  }
  for(l in 1:(p - 1)){
    for(m in (l + 1):p){
      # grid.distance = sqrt(sum((coor[l, ] - coor[m, ]) ^ 2))
      grid.distance = min(coor[j] - coor[i], p - coor[j] + coor[i])
      if(grid.distance <= bw){
        Sigma2_lm = sum(x.scale[l, ] * x.scale[m, ])
        for(i in 1:n){
          mapping.matrix[i, i] = mapping.matrix[i, i] + 2 * x.scale.inv[i, l] * x.scale.inv[i, m] * Sigma2_lm
        }
        for(i in 1:(n - 1)){
          for(j in (i + 1):n){
            mapping.matrix[i, j] = mapping.matrix[j, i] = mapping.matrix[i, j] + (x.scale.inv[i, l] * x.scale.inv[j, m] + x.scale.inv[i, m] * x.scale.inv[j, l]) * Sigma2_lm
          }
        }
      }
    }
  }
  
  mapping.matrix.eigen = eigen(mapping.matrix)
  mapping.matrix.eigen.value = mapping.matrix.eigen$values
  mapping.matrix.eigen.vector = mapping.matrix.eigen$vectors
  mapping.matrix.eigen.value[which(mapping.matrix.eigen.value < 1e-5)] = 1e-5
  z = x.scale %*% mapping.matrix.eigen.vector %*% diag(sqrt(mapping.matrix.eigen.value))
  
  return(z)
}