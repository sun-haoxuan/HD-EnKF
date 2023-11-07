bandP = function(K, N, is.circ = T){
  Bk = diag(1, N)
  if(K <= 0){
    Bk = diag(1, N)
  }else if(K < N){
    idx1 = 1:N
    for(k in 1:K){
      if(is.circ){
        idx2 = (idx1 + k) %% N
        idx2[which(idx2 == 0)] = N
        for(j in 1:N){
          Bk[idx1[j], idx2[j]] = Bk[idx2[j], idx1[j]] = 1
        }
      }else{
        idx2 = idx1[1:(N - k)] + k
        for(j in 1:(N - k)){
          Bk[idx1[j], idx2[j]] = Bk[idx2[j], idx1[j]] = 1
        }
      }
    }
  }else{
    for(i in 1:N){
      for(j in 1:N){
        Bk[i, j] = 1
      }
    }
  }
  return(Bk)
}

taperP = function(k, N, is.circ = T){
  Bk = matrix(0, N, N)
  if(k <= 0){
    Bk = diag(1, N)
  }else if(k < N){
    for(i in 1:N){
      for(j in i:N){
        if(is.circ){
          l = min(j - i, N - (j - i))
        }else{
          l = j - i
        }
        if(l <= k){
          Bk[i, j] = Bk[j, i] = 1
        }else if(k < l & l < 2 * k){
          Bk[i, j] = Bk[j, i] = 2 - l / k
        }
      }
    }
  }else{
    for(i in 1:N){
      for(j in 1:N){
        Bk[i, j] = 1
      }
    }
  }
  return(Bk)
}


