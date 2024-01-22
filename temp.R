setwd("E:/Project/EnKF_GTE")
source('Code/Functions/cal_taper_cov_L96.R')
plot1 = function(x){corrplot::corrplot(x, method = 'color', is.corr = F)}

Z = (x.f - x.f.bar) / sqrt(n - 1)
bw = cal_taper_bandwidth(t(Z), mat.dist, BL, interval)
x = t(Z)

n = nrow(x)
p = ncol(x)

Sn = t(x) %*% x
mat.dist.taper = diag(1, p)
for(i in 1:(p - 1)){
  for(j in (i + 1):p){
    mat.dist.taper[i, j] = mat.dist.taper[j, i] = BL(mat.dist[i, j] / bw)
  }
}

Sigma.hat = Sn * mat.dist.taper

Sigma.hat.eigen = RSpectra::eigs(Sigma.hat, p)
Sigma.hat.eigen.values = Sigma.hat.eigen$values
Sigma.hat.eigen.vectors = Sigma.hat.eigen$vectors
length(which(Sigma.hat.eigen.values > 1e-5))
Sigma.hat.eigen.values[which(Sigma.hat.eigen.values < 1e-5)] = 1e-5
Sigma.tilde = Sigma.hat.eigen.vectors %*% diag(Sigma.hat.eigen.values) %*% t(Sigma.hat.eigen.vectors)


mat.dist.taper.eigen = RSpectra::eigs(mat.dist.taper, p)
mat.dist.taper.eigen.values = mat.dist.taper.eigen$values
mat.dist.taper.eigen.vectors = mat.dist.taper.eigen$vectors
mat.dist.taper.eigen.values[which(mat.dist.taper.eigen.values < 1e-5)] = 1e-5
mat.dist.taper.tilde = mat.dist.taper.eigen.vectors %*% diag(mat.dist.taper.eigen.values) %*% t(mat.dist.taper.eigen.vectors)
Sigma.module = Sn * mat.dist.taper.tilde
Sigma.module.eigen = RSpectra::eigs(Sigma.module, p)

plot(Sigma.hat.eigen.values, Sigma.module.eigen$values)

plot1(Sigma.hat)
plot1(Sigma.tilde)
plot1(Sigma.module)
plot1(Sigma.tilde - Sigma.module)
