cal_bandwidth_L96_old = function(x, mu = NULL){
  p = nrow(x)
  n = ncol(x)
  if(is.null(mu)){
    x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  }else{
    x.scale = (x - mu) / sqrt(n - 1)
  }
  # coor = matrix(c(cos(1:p * 2 * pi / p), sin(1:p * 2 * pi / p)), p, 2)
  coor = 1:p
  
  M_n = function(k){
    tmp = 0
    for(i in 1:(p - 1)){
      sigma_ii = sum(x.scale[i, ] ^ 2)
      for(j in i:p){
        # grid.distance = sqrt(sum((coor[j, ] - coor[i, ]) ^ 2))
        grid.distance = min(coor[j] - coor[i], p - coor[j] + coor[i])
        if(grid.distance > k){
          tmp = tmp + sum(x.scale[i, ] * x.scale[j, ]) / p
        }else{
          tmp = tmp + sigma_ii * sum(x.scale[j, ] ^ 2) / p / n
        }
      }
    }
    return(tmp)
  }
  
  opt = optimise(M_n, c(0, p / 2), tol = 1)
  bw = opt$minimum
  return(bw)
}