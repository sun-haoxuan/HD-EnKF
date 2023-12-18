banding = function(z){
  if(z <= 1){
    return(1)
  }else{
    return(0)
  }
}

tapering = function(z){
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

p = 100
n = 30

matrix.distance = diag(0, p)
nu = 10
nv = p / 10
grid = data.frame(u = rep(1:nu-1, each = nv), v =rep(1:nv-1, times = nu))
for(i in 1:(p - 1)){
  for(j in (i + 1):p){
    matrix.distance[i, j] = matrix.distance[j, i] = round(sqrt(sum((grid[i, ] - grid[j, ]) ^ 2)), 1)
  }
}
Sigma = 0.5 ^ matrix.distance
eg = eigen(Sigma)
eg.value = eg$values
eg.value[which(eg.value < 1e-5)] = 1e-5
Sigma.sqrt = eg$vectors %*% diag(eg.value ^ 0.5) %*% t(eg$vectors)

# corrplot::corrplot(Sigma.GTE.30, method = 'color', is.corr = F)
# data.frame(
#   x = p - rep(1:p, times = p),
#   y = rep(1:p, each = p),
#   col = as.vector(Sigma)
# ) %>% ggplot()+
#   geom_tile(aes(x = x, y = y, fill = col))+
#   scale_fill_gradient(low = '#FFFFFF',  high = '#FF0000', limits = c(0, 1))+
#   theme_bw()+
#   theme(
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     # legend.key.height = unit(1, 'in'),
#     legend.title = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_blank()
#   )

result = foreach(mt = c('banding', 'tapering', 'GC'), .combine = 'rbind') %do%
  {
    if(mt == 'banding'){taper.func = banding}
    if(mt == 'tapering'){taper.func = tapering}
    if(mt == 'GC'){taper.func = GC}
    
    set.seed(1234)
    x = Sigma.sqrt %*% matrix(rnorm(p * n, 0, 1), p, n)
    
    dist.all = sort(unique(as.vector(matrix.distance)))
    risk = list()
    for(d_ij in dist.all){
      list.grid = which(matrix.distance == d_ij, arr.ind = TRUE)
      tmp1 = tmp2 = 0
      for(grd in 1:nrow(list.grid)){
        i = list.grid[grd, 1]
        j = list.grid[grd, 2]
        x_i_2 = x[i, ] ^ 2
        x_j_2 = x[j, ] ^ 2
        x_i_j = x[i, ] * x[j, ]
        xi2_xj2 = sum(x_i_2 * x_j_2)
        tmp1 = tmp1 + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
        tmp2 = tmp2 + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2)
      }
      risk[[as.character(d_ij)]] = c(tmp1, tmp2)
    }
    risk_taper = function(k){
      list.dist = names(risk)[order(as.numeric(names(risk)))]
      tmp = 0
      id_dist = max(which(as.numeric(list.dist) <= k))
      if(!is.infinite(id_dist)){
        for(i in 1:id_dist){
          # g_ij = banding(as.numeric(list.dist)[i] / k)
          # g_ij = tapering(as.numeric(list.dist)[i] / k)
          # g_ij = GC(as.numeric(list.dist)[i] / k)
          g_ij = taper.func(as.numeric(list.dist)[i] / k)
          tmp = tmp + (g_ij ^ 2 - 2 * g_ij) * risk[[list.dist[i]]][1] + g_ij ^ 2 * risk[[list.dist[i]]][2] / n
        }
      }
      return(tmp)
    }
    opt.GTE = optimise(risk_taper, c(0, max(matrix.distance)))
    bw.GTE = opt.GTE$minimum
    
    Sigma.sample = cov(t(x))
    mat.dist.taper = 0.5*diag(1, p)
    for(i in 1:(p - 1)){
      for(j in (i + 1):p){
        mat.dist.taper[i, j] = mat.dist.taper[j, i] = taper.func(matrix.distance[i, j] / bw.GTE)
      }
    }
    Sigma.GTE = mat.dist.taper * cov(t(x))
    
    eg = eigen(Sigma.GTE)
    eg.value = eg$values
    eg.value[which(eg.value < 1e-5)] = 1e-5
    Sigma.GTE.10 = eg$vectors[, 1:10] %*% diag(eg.value[1:10]) %*% t(eg$vectors[, 1:10])
    Sigma.GTE.30 = eg$vectors[, 1:30] %*% diag(eg.value[1:30]) %*% t(eg$vectors[, 1:30])
    Sigma.GTE.50 = eg$vectors[, 1:50] %*% diag(eg.value[1:50]) %*% t(eg$vectors[, 1:50])
    Sigma.GTE.60 = eg$vectors[, 1:60] %*% diag(eg.value[1:60]) %*% t(eg$vectors[, 1:60])
    Sigma.GTE.70 = eg$vectors[, 1:70] %*% diag(eg.value[1:70]) %*% t(eg$vectors[, 1:70])
    Sigma.GTE.90 = eg$vectors[, 1:90] %*% diag(eg.value[1:90]) %*% t(eg$vectors[, 1:90])
    Sigma.GTE.100 = eg$vectors[, 1:100] %*% diag(eg.value[1:100]) %*% t(eg$vectors[, 1:100])
    
    
    data.frame(
      method = mt,
      bw = bw.GTE,
      type = c('Sn', 'origin', 'ev10', 'ev30', 'ev50', 'ev60', 'ev70', 'ev90', 'ev100'),
      rawdiff = c(norm(Sigma - Sigma.sample, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE.10, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE.30, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE.50, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE.60, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE.70, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE.90, 'F') / sqrt(p),
                  norm(Sigma - Sigma.GTE.100, 'F') / sqrt(p)),
      invdiff = c(norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.sample + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE.10 + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE.30 + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE.50 + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE.60 + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE.70 + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE.90 + 0.5*diag(p)), 'F') / sqrt(p),
                  norm(solve(Sigma + 0.5*diag(p))- solve(Sigma.GTE.100 + 0.5*diag(p)), 'F') / sqrt(p))
    )
  }
result$type = factor(result$type, c('Sn', 'origin', 'ev10', 'ev30', 'ev50', 'ev60', 'ev70', 'ev90', 'ev100'))
result$method = factor(result$method, c('banding', 'tapering', 'GC'))
result %>%
  filter(type %in% c('Sn', 'ev10', 'ev30', 'ev50', 'ev70', 'ev90')) %>%
  dplyr::select(method, type, rawdiff) %>%
  spread(type, rawdiff) %>% 
  xtable::xtable(digits = 3) %>% print(include.rownames = F)
result %>%
  filter(type %in% c('Sn', 'ev10', 'ev30', 'ev50', 'ev70', 'ev90')) %>%
  dplyr::select(method, type, invdiff) %>%
  spread(type, invdiff) %>% 
  xtable::xtable(digits = 3) %>% print(include.rownames = F)

