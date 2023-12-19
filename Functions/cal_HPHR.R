HPHR.det = function(lambda, D, R){
  return(det(R) * prod(lambda * D ^ 2 + 1))
}

HPHR.det2 = function(lambda, D, R, shrink){
  return(det(R) * prod((lambda * D ^ 2 + 1) / shrink))
}

HPHR.inv = function(lambda, HZ, R.inv, n){
  tmp = solve(diag(n) + lambda * t(HZ) %*% R.inv %*% HZ) %*% t(HZ) %*% R.inv
  return(R.inv - lambda * R.inv %*% HZ %*% tmp)
}
