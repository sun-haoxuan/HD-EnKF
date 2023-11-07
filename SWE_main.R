rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)
library(Rcpp)
library(microbenchmark)
library(RSpectra)

setwd("E:/Project/EnKFTapering")
dir.open = 'Code/'
dir.save = 'Output/SWE/SWE_k5e4_sigma0.01/'

source(paste0(dir.open, 'Functions/cal_HPHR.R'))
source(paste0(dir.open, 'Functions/cal_RMSE.R'))
source(paste0(dir.open, 'Functions/SWE_Forecast.R'))
source(paste0(dir.open, 'Functions/SWE_State.R'))
source(paste0(dir.open, 'Functions/SWE_Analyse.R'))

param = list()
for(mt in c('freerun', 'oracle', 'inflation', 'banding', 'tapering', 'threshold')){
  for(k in c(1e4, 5e4)){
    for(b in 1:1){
      if(mt %in% c('freerun', 'oracle', 'inflation')){
        ea = iu = F
      }else{
        ea = iu = T
      }
      param[['k']] = c(param[['k']], k)
      param[['b']] = c(param[['b']], b)
      param[['mt']] = c(param[['mt']], mt)
      param[['ea']] = c(param[['ea']], ea)
      param[['iu']] = c(param[['iu']], iu)
    }
  }
}
param = as.data.frame(param)

dir.create(dir.save)
dir.create(paste0(dir.save, 'temp'))
for(mt in unique(param$mt)){
  dir.create(paste0(dir.save, mt))
  for(k in unique(param$k)){
    dir.create(paste0(dir.save, mt, '/K', k))
  }
}

cores = 10
cl = makeCluster(cores)
registerDoParallel(cl, cores = cores)

t1 = Sys.time()
result = foreach(
  i = 1:nrow(param),
  .packages = c('tidyverse', 'microbenchmark', 'Rcpp', 'RSpectra'),
  .combine = 'rbind') %do%
  {
    print(param[i, ])
    k = param$k[i]
    b = param$b[i]
    mt = param$mt[i]
    ea = param$ea[i]
    iu = param$iu[i]
    
    set.seed(1234 + b)
    opt = list(
      ## Generate true state
      nx = 50, # grid size of x
      ny = 31, # grid size of ys
      L = 5e+5,
      D = 3e+5,
      g = 9.8, # gravitational constant
      f = 1e-4,
      k = 1e+4,
      dt = 30, # hardwired timestep
      h0 = 50,
      h1 = 5.5,
      h2 = 3.325,
      ## Generate observations
      id = sort(sample(50, 10)),
      sigma.obs = c(sqrt(0.5), sqrt(0.5), 1), # observation error
      rho = 0.5, ## observation
      vf = 1, # for Wrong R
      ## Assimilation Step
      S.o.seq = 360, # assimilate every 3h
      S.o.start = 61, ## Start of the observation
      S.o.end = 5820, # assimilate 48h, then forecast 24h
      S = 5820,
      ## Assimilation
      m = 100,
      k.set = k,
      sigma.i = 1,
      method = mt,
      eig_adj = ea,
      iter_update = iu
    )
    
    if(mt == 'oracle'){
      opt$m = 1000
    }
    
    state = SWE_State(opt)
    analyse = tryCatch(
      expr = {
        SWE_Analyse(state, opt, full_output = F)
      },
      error = function(e){
        tryCatch(
          expr = {
            SWE_Analyse(state, opt, full_output = F)
          },
          error = function(e){
            list(converge = 1, error = rep(NA, 5820))
          }
        )
      }
    )
    save(analyse, file = paste0(dir.save, mt, '/K', k, '/B', b, '_n1000.Rdata'))

    # load(paste0(dir.save, mt, '/K', k, '/B', b, '.Rdata'))
    
    data.frame(
      method = mt,
      k.set = k,
      boot = b,
      cv = analyse$converge,
      RMSE = sqrt(mean(analyse$error[2941:5820] ^ 2))
    )
  }
stopCluster(cl)
stopImplicitCluster()

result %>%
  group_by(k.set, method) %>%
  summarise(ave = mean(RMSE, na.rm = T)) %>%
  spread(method, ave) %>%
  print()

t2 = Sys.time()
print(t2 - t1)

