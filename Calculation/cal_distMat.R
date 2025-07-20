cal_distMat = function(option){
  p = option$p

  distance.matrix = diag(0, p)
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      distance.matrix[i, j] = distance.matrix[j, i] = min(j - i, p - j + i)
    }
  }
  
  return(distance.matrix)
}