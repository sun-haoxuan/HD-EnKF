cal_bandwidth_L96 = function(x){
  p = nrow(x)
  n = ncol(x)
  x.scale = (x - apply(x, 1, mean)) / sqrt(n - 1)
  coor = matrix(c(cos(1:p * 2 * pi / p), sin(1:p * 2 * pi / p)), p, 2)
  
  M_n = function(k){
    tmp = 0
    for(i in 1:(p - 1)){
      sigma_ii = sqrt(sum(x.scale[i, ] ^ 2))
      for(j in i:p){
        grid.distance = sqrt(sum((coor[j, ] - coor[i, ]) ^ 2))
        if(grid.distance > k){
          tmp = tmp + sum(x.scale[i, ] * x.scale[j, ]) / p
        }else{
          tmp = tmp + sigma_ii * sqrt(sum(x.scale[j, ] ^ 2)) / p / n
        }
      }
    }
    return(tmp)
  }
  
  opt = optimise(M_n, c(0, 2), tol = 1)
  bw = opt$minimum
  return(bw)
}