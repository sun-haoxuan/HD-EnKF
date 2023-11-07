cal_GTE_bandwidth_SWE_version1 = function(x, dist.func, taper.func, digits = 1, limit = Inf, tol = .Machine$double.eps^0.25){
  p = ncol(x)
  n = nrow(x)
  coor = data.frame(
    x = rep(rep(1:50, times = 31), 3),
    y = rep(rep(1:31, each = 50), 3)
  )
  
  # t1 = Sys.time()
  risk = list()
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      # d_ij = round(dist.func(i, j, coor), digits)
      d_ij = round(sqrt(sum((coor[i, ] - coor[j, ]) ^ 2)), digits)
      if(d_ij <= limit){
        d_ij = as.character(d_ij)
        if(! d_ij %in% names(risk)){
          risk[[d_ij]] = c(0, 0)
        }
        x_i_2 = x[, i] ^ 2
        x_j_2 = x[, j] ^ 2
        x_i_j = x[, i] * x[, j]
        xi2_xj2 = sum(x_i_2 * x_j_2)
        risk[[d_ij]][1] = risk[[d_ij]][1] + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
        risk[[d_ij]][2] = risk[[d_ij]][2] + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
      }
    }
  }
  # t2 = Sys.time()
  # print(t2 - t1)
  # save(risk, file = 'test_risk_oldcode2.Rdata')
  
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
  
  opt = optimise(function(k){risk_taper}, c(0, max(as.numeric(names(risk)))), tol = tol)
  
  return(opt$minimum)
}

cal_GTE_bandwidth_SWE = function(x, dist.func, taper.func, digits = 1, limit = Inf, tol = .Machine$double.eps^0.25){
  p = ncol(x)
  n = nrow(x)
  coor = data.frame(
    x = rep(rep(1:50, times = 31), 3),
    y = rep(rep(1:31, each = 50), 3)
  )
  
  risk = list()
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      # d_ij = round(dist.func(i, j, coor), digits)
      d_ij = round(sqrt(sum((coor[i, ] - coor[j, ]) ^ 2)), digits)
      if(d_ij <= limit){
        d_ij = as.character(d_ij)
        if(! d_ij %in% names(risk)){
          risk[[d_ij]] = c(0, 0)
        }
        x_i_2 = x[, i] ^ 2
        x_j_2 = x[, j] ^ 2
        x_i_j = x[, i] * x[, j]
        xi2_xj2 = sum(x_i_2 * x_j_2)
        risk[[d_ij]][1] = risk[[d_ij]][1] + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
        risk[[d_ij]][2] = risk[[d_ij]][2] + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
      }
    }
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
  
  opt = optimise(function(k){risk_taper}, c(0, max(as.numeric(names(risk)))), tol = tol)
  
  return(opt$minimum)
}