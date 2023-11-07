ev_adjust = function(P, m){
  require(RSpectra)
  
  ev = eigs(P, m, which = "LM")
  A1 = ev$vectors
  A2 = ev$values
  A2[which(A2 < 1e-5)] = 0
  P.adj = A1 %*% diag(A2) %*% t(A1)
  
  return(P.adj)
}
