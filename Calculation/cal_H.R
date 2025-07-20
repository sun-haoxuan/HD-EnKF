cal_H = function(option){
  set.seed(option$boot)
  
  p = option$p
  q = option$q
  
  p.o = sort(sample(p, q))
  nsite = length(p.o)
  
  H = matrix(0, nsite, p)
  for(k in 1:nsite){
    H[k, p.o[k]] = 1
  }
  
  return(H)
}