bw =cal_taper_bandwidth(t(z), mat.dist, banding, interval)

z1 = cal_taper_Z(t(z), bw, mat.dist, banding)
P.hat.1 = z1 %*% t(z1)



Sigma.hat = Z %*% t(Z)
mat.dist.taper = diag(1, p)
for(i in 1:(p - 1)){
  for(j in (i + 1):p){
    mat.dist.taper[i, j] = mat.dist.taper[j, i] = banding(mat.dist[i, j] / bw)
  }
}
Sigma.hat = Sigma.hat * mat.dist.taper

A = svd(z)
U = A$u
d = A$d
V = A$v
# Z - U %*% diag(d) %*% t(V)
Z.inv = V[, 1:(n-1)] %*% diag(d[1:(n-1)] ^(-1)) %*% t(U[, 1:(n-1)])
AAT = Z.inv %*% Sigma.hat %*% t(Z.inv)
A1 = eigen(AAT)
eg.vec = A1$vectors
eg.value = A1$values
AAT - eg.vec %*% diag(eg.value) %*% t(eg.vec)

eg.value[which(eg.value < 1e-5)] = 1e-5
A.hat = eg.vec %*% diag(eg.value ^ (0.5)) %*% t(eg.vec)
AAT - A.hat %*% t(A.hat)

Z.tilde = Z %*% A.hat
P.hat.2 = Z.tilde %*% t(Z.tilde)

AAT.hat = matrix(NA, n, n)
for(j in 1:n){
  for(k in j:n){
    tmp = 0
    for(l1 in 1:p){
      for(l2 in 1:p){
        tmp = tmp + Z.inv[j, l1] * Z.inv[k, l2] * sum(Z[l1, ] * Z[l2, ]) * mat.dist.taper[l1, l2]
      }
    }
    AAT.hat[j, k] = AAT.hat[k, j] = tmp
  }
}

corrplot::corrplot(t(x.tilde) %*% x.tilde, method = 'color', is.corr = F)
