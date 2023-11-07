rm(list = ls())
setwd("E:/Project/EnKF")
library(tidyverse)
library(foreach)
library(doParallel)
library(MASS)
library(Rcpp)
library(microbenchmark)
library(RSpectra)

source('Code/Functions/band_Mat.R')
source('Code/Functions/cal_HPHR.R')
source('Code/Functions/cal_RMSE.R')
source('Code/Functions/cal_threshold.R')
source('Code/Functions/eig_adj.R')
source('Code/Functions/L96_State.R')
source('Code/Functions/L96_Analyse.R')

param = list()
# for(mt in c('freerun', 'origin', 'oracle', 'inflation', 'banding', 'tapering', 'threshold')){
for(mt in c('banding', 'tapering')){
  for(f in 4:12){
    for(b in 1:100){
      if(mt %in% c('freerun', 'origin', 'oracle', 'inflation')){
        ea = iu = F
      }else{ 
        ea = iu = T
      }
      param[['f']] = c(param[['f']], f)
      param[['b']] = c(param[['b']], b)
      param[['mt']] = c(param[['mt']], mt)
      param[['ea']] = c(param[['ea']], ea)
      param[['iu']] = c(param[['iu']], iu)
    }
  }
}
param = as.data.frame(param)

setting = data.frame(
  p = NA,
  q = NA,
  n = NA,
  sigma.tr = NA
)
setting = setting %>%
  # rbind(c(40, 40, 30, 0)) %>%
  # rbind(c(40, 20, 30, 0)) %>%
  # rbind(c(40, 30, 30, 0)) %>%
  # rbind(c(100, 100, 30, 0)) %>%
  # rbind(c(100, 100, 60, 0)) %>%
  rbind(c(100, 100, 90, 0)) %>%
  slice(-1)

cores = 10
cl = makeCluster(cores)
registerDoParallel(cl, cores = cores)

for(i in 1:nrow(setting)){
  print(setting[i, ])
  p = setting$p[i]
  q = setting$q[i]
  n = setting$n[i]
  sigma.tr = setting$sigma.tr[i]
  
  dir.save = paste0('Output/L96/Data_1001/L96_p', p, '_q', q, '_n', n, '_sigma', sigma.tr, '/')
  dir.create(dir.save)
  for(mt in unique(param$mt)){
    dir.create(paste0(dir.save, mt))
    for(f in unique(param$f)){
      dir.create(paste0(dir.save, mt, '/F', f))
    }
  }
  
  t1 = Sys.time()
  result = foreach(
    j = 1:nrow(param),
    .packages = c('tidyverse', 'MASS', 'Rcpp', 'microbenchmark', 'RSpectra'),
    .combine = 'rbind'
  ) %dopar% {
    sourceCpp('Code/Functions/cal_bandwidth_m0.cpp')
    
    f = param$f[j]
    b = param$b[j]
    mt = param$mt[j]
    ea = param$ea[j]
    iu = param$iu[j]
    filename = paste0(mt, '_F', f, '_B', b, '.Rdata')
    
    set.seed(1234 + b)
    
    opt = list(
      p = p,
      S = 2000,
      F.true = 8,
      h = 0.05,
      sigma.tr = sigma.tr,
      p.o = sort(sample(p, q)),
      S.o = seq(4, 2000, 4),
      sigma = 1,
      rho = 0.5,
      n = n,
      sigma.in = 0.1,
      F.set = f,
      method = mt,
      eig_adj = ea,
      iter_update = iu,
      is.circ = T
    )
    
    if(mt == 'oracle'){
      opt$n = 1000
    }
    if(q == 20){
      opt$S.o = seq(2, 2000, 2)
    }
    if(q == 30){
      opt$S.o = seq(3, 2000, 3)
    }
    
    state = L96_State(opt)
    analyse = tryCatch(
      expr = {
        L96_Analyse(state, opt, full_output = F)
      },
      error = function(e){
        tryCatch(
          expr = {
            L96_Analyse(state, opt, full_output = F)
          },
          error = function(e){
            list(converge = 1, error = rep(NA, 2000))
          }
        )
      }
    )
    save(analyse, file = paste0(dir.save, mt, '/F', f, '/', filename))
    
    # load(paste0(dir.save, mt, '/F', f, '_B', b, '.Rdata'))
    
    data.frame(
      method = mt,
      F.set = f,
      boot = b,
      cv = analyse$converge,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2)),
      obj = mean(analyse$objective[251:500])
    )
  }

  result %>%
    group_by(F.set, method) %>%
    summarise(ave = mean(RMSE, na.rm = T)) %>%
    spread(method, ave) %>%
    print()
  
  t2 = Sys.time()
  print(t2 - t1)
}
stopCluster(cl)
stopImplicitCluster()
