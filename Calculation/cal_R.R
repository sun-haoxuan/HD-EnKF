cal_R = function(option){
  set.seed(option$boot)
  
  p = option$p
  q = option$q
  sigma = option$sigma
  rho = option$rho
  
  p.o = sort(sample(p, q))
  nsite = length(p.o)

  R = matrix(0, nsite, nsite)
  for(i in 1:nsite){
    for(j in 1:nsite){
      R[i,j] = sigma ^ 2 * rho ^ min(abs(i - j), nsite - abs(i - j))
    }
  }
  R.inv = solve(R)
  R.ev = eigen(R)
  R.inv.root = R.ev$vectors %*% diag(R.ev$values ^ (-0.5)) %*% t(R.ev$vectors)
  
  output = list(
    R = R,
    R.inv = R.inv,
    R.inv.root = R.inv.root
  )
  
  return(output)
}