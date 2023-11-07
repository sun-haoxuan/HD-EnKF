rm(list = ls())
setwd("E:/Project/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)

source('Code/Functions/cal_HPHR.R')
source('Code/Functions/cal_RMSE.R')
source('Code/Functions/SWE_State.R')
source('Code/Functions/SWE_Forecast.R')
source('Code/Functions/SWE_Analyse.R')
source('Code/cal_taper_cov_SWE.R')

cores = 10
cl = makeCluster(cores, outfile = 'Output/Result_SWE_1103/outfile.txt')
registerDoParallel(cl, cores = cores)

result = foreach(b = 1:10, .combine = 'rbind') %dopar%
  {
    opt = list(
      ## Generate true state
      nx = 50, # grid size of x
      ny = 31, # grid size of ys
      L = 5e+5,
      D = 3e+5,
      g = 9.8, # gravitational constant
      f = 1e-4,
      k = NA,
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
      S = 8700,
      ## Assimilation
      n = 100,
      k.true = 1e+4,
      k.set = 5e+4,
      sigma.i = 1,
      method = 'inflation',
      boot = b
    )
    
    set.seed(1234 + b)
    state = SWE_State(opt)
    analyse = tryCatch(
      expr = {
        SWE_Analyse(state, opt)
      },
      error = function(e){
        list(converge = 1, error = matrix(NA, 4, 8700))
      }
    )
    
    save(analyse, file = paste0('Output/Result_SWE_1103/Inflation_B', b, '.Rdata'))
    
    data.frame(
      method = 'inflation',
      boot = b,
      cv = analyse$converge,
      RMSE_pre = sqrt(mean(analyse$error[1, 2941:5820] ^ 2)),
      RMSE_post = sqrt(mean(analyse$error[1, 5821:8700] ^ 2)),
      RMSE_u_pre = sqrt(mean(analyse$error[2, 2941:5820] ^ 2)),
      RMSE_u_post = sqrt(mean(analyse$error[2, 5821:8700] ^ 2)),
      RMSE_v_pre = sqrt(mean(analyse$error[3, 2941:5820] ^ 2)),
      RMSE_v_post = sqrt(mean(analyse$error[3, 5821:8700] ^ 2)),
      RMSE_h_pre = sqrt(mean(analyse$error[4, 2941:5820] ^ 2)),
      RMSE_h_post = sqrt(mean(analyse$error[4, 5821:8700] ^ 2))
    )
  }

stopCluster(cl)
stopImplicitCluster()
