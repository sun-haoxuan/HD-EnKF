rm(list = ls())
setwd("E:/Project/EnKF")
library(tidyverse)
library(foreach)
library(Rcpp)

sourceCpp('Code/Functions/cal_bandwidth_m0.cpp')

set.p = c(40, 100, 500, 1000)
set.n = c(30, 100)
result = foreach(
  p = rep(set.p, each = length(set.n)),
  n = rep(set.n, times = length(set.p)),
  .combine = 'rbind'
  ) %do%
  {
    print(paste(p, n))
    theta = 1 / 0.9
    Sigma = matrix(0, p, p)
    for(l1 in 1:p){
      for(l2 in 1:p){
        Sigma[l1, l2] = theta ^ (-abs(l1 - l2))
      }
    }
    eg = eigen(Sigma)
    Sigma.sqrt = eg$vectors %*% diag(eg$values ^ 0.5) %*% t(eg$vectors)
    x = Sigma.sqrt %*% matrix(rnorm(p * n, 0, 1), p, n)
    
    t1 = Sys.time()
    bw = Ccal_banding(x, p - 1)
    t2 = Sys.time()
    
    data.frame(
      p = p,
      n = n,
      bw = bw,
      t1 = t1,
      t2 = t2,
      time = t2 - t1
    )
  }
write.csv(result, 'result_test_bandwidth.csv')
