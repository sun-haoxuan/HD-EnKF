rm(list = ls())
setwd("E:/Project/EnKFTapering")

library(tidyverse)
library(foreach)

source('Code/Functions/cal_GTE_bandwidth.R')
source('Code/Functions/dist_functions.R')
source('Code/Functions/taper_functions.R')
Rcpp::sourceCpp('Code/Functions/cal_bandwidth_m0.cpp')

## Generate data ####
theta = 0.9
p = 1000
n = 100
Sigma = matrix(NA, p, p)
for(i in 1:p){
  for(j in 1:p){
    Sigma[i, j] = theta ^ (abs(i - j))
  }
}
Sigma.eigen = eigen(Sigma)
Sigma.root = Sigma.eigen$vectors %*% diag(Sigma.eigen$values ^ 0.5) %*% t(Sigma.eigen$vectors)
z = matrix(rnorm(p * n), p, n)
x.raw = Sigma.root %*% z
x = t(x.raw - apply(x.raw, 1, mean))

# Ccal_banding(x, n - 1)
Cband.loss = Cbandloss(x, n - 1)
which.min(Cband.loss)
plot(Cband.loss, type = 'l')

# cal_GTE_bandwidth(x, dist_dim1, banding, tol = 1)
risk = risk_each_dist(x, dist_dim1)
prop.loss = foreach(k = 0:(n - 1), .combine = c) %do% {risk_taper(risk, k, banding)}
which.min(prop.loss)

dist.func = dist_dim1
taper.func = banding

t1 = Sys.time()
tmp = list()
for(i in 1:(p - 1)){
  for(j in (i + 1):p){
    # print(paste(i, j))
    d_ij = dist.func(i, j, p)
    if(! d_ij %in% names(tmp)){
      tmp[[d_ij]] = c(0, 0)
    }
    # for(l in 1:n){
    #   for(m in setdiff(1:n, l)){
    #     tmp[[d_ij]][1] = tmp[[d_ij]][1] + (x[l, i] ^ 2) * (x[m, j] ^ 2) / n
    #     # tmp[[d_ij]][2] = tmp[[d_ij]][2] + x[l, i] * x[l, j] * x[m, i] * x[m, j]
    #   }
    # }
    x_i_2 = x[, i] ^ 2
    x_j_2 = x[, j] ^ 2
    x_i_j = x[, i] * x[, j]
    xi2_xj2 = sum(x_i_2 * x_j_2)
    tmp[[d_ij]][1] = tmp[[d_ij]][1] + (sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
    tmp[[d_ij]][2] = tmp[[d_ij]][2] + sum(x_i_j %*% t(x_i_j)) - xi2_xj2
  }
}
t2 = Sys.time()
print(t2 - t1)

# foreach(i = 1:199, .combine = 'rbind') %do% {data.frame(true = tmp0[[i]][1], new = tmp[[i]][1])} %>% plot()

i = 1; j = 2; a = b = 0
for(l in 1:n){
  for(m in setdiff(1:n, l)){
    a = a + (x[l, i] ^ 2) * (x[m, j] ^ 2) / n
    b = b + x[l, i] * x[l, j] * x[m, i] * x[m, j]
  }
}

x_i_2 = x[, i] ^ 2
x_j_2 = x[, j] ^ 2
x_i_j = x[, i] * x[, j]
xi2_xj2 = sum(x_i_2 * x_j_2)
(sum(x_i_2 %*% t(x_j_2))- xi2_xj2) / n
sum(x_i_j %*% t(x_i_j)) - xi2_xj2



