dist_dim1 = function(i, j, coor){return(j - i)}

dist_dim1_circ = function(i, j, p){return(min(j - i, p - j + i))}

## SWE ####
coor = data.frame(
  x = rep(rep(1:50, times = 31), 3),
  y = rep(rep(1:31, each = 50), 3)
)
dist_dim2 = function(i, j, coor){
  return(sqrt(sum((coor[i, ] - coor[j, ]) ^ 2)))
}
risk_each_dist = function(x, dist.func, digits = 0, limit = Inf){
  p = ncol(x)
  n = nrow(x)
  
  digits = 0
  limit = 30
  t1 = Sys.time()
  tmp = list()
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      # print(paste(i, j))
      d_ij = round(dist.func(i, j, coor), digits)
      if(d_ij <= limit){
        d_ij = as.character(d_ij)
        if(! d_ij %in% names(tmp)){
          tmp[[d_ij]] = c(0, 0)
        }
        x_i_2 = x[, i] ^ 2
        x_j_2 = x[, j] ^ 2
        x_i_j = x[, i] * x[, j]
        xi2_xj2 = sum(x_i_2 * x_j_2)
        tmp[[d_ij]][1] = tmp[[d_ij]][1] + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
        tmp[[d_ij]][2] = tmp[[d_ij]][2] + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
      }
    }
  }
  risk1 = tmp
  t2 = Sys.time()
  print(paste(digits, limit, t2 - t1))
  return(tmp)
}
