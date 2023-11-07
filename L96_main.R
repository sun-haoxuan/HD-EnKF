rm(list = ls())
setwd("/disk/home/sunhaoxuan/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)

source('Code/Functions/cal_GTE_bandwidth.R')
source('Code/Functions/cal_HPHR.R')
source('Code/Functions/cal_RMSE.R')
source('Code/Functions/ev_adjust.R')
source('Code/Functions/L96_State.R')
source('Code/Functions/L96_Analyse.R')
source('Code/Functions/RK4.R')
source('Code/Functions/taper_functions.R')
source('Code/Functions/dist_functions.R')
source('Code/Functions/cal_taper_cov_L96.R')

methods = c('oracle', 'inflation & iteration', 'GC', 'inflation then GC')
param = list()
for(mt in methods){
  for(f in 4:12){
    for(b in 1:1){
      param[['f']] = c(param[['f']], f)
      param[['b']] = c(param[['b']], b)
      param[['mt']] = c(param[['mt']], mt)
    }
  }
}
param = as.data.frame(param)

dir.save = 'Output/L96_1025_p300/'
dir.create(dir.save)
for(mt in unique(param$mt)){
  dir.create(paste0(dir.save, mt))
}

cores = 36
cl = makeCluster(cores)
registerDoParallel(cl, cores = cores)
result = foreach(j = 1:nrow(param),  .packages = c('Matrix', 'RSpectra'),.combine = 'rbind') %dopar%
  {
    print(param[j, ])
    f = param$f[j]
    b = param$b[j]
    mt = param$mt[j]
    filename = paste0(mt, '_F', f, '_B', b, '.Rdata')
    
    set.seed(1234 + b)
    option = list(
      p = 300,
      S = 2000,
      F.true = 8,
      h = 0.05,
      p.o = seq(1, 300, 1),
      S.o = seq(4, 2000, 4),
      sigma = 1,
      rho = 0.5,
      n = 30,
      method = mt,
      F.set = f,
      sigma.in = f / 40
    )
    state = L96_State(option)
    analyse = tryCatch(
      expr = {
        L96_Analyse(state, option)
      },
      error = function(e){
        list(converge = 1, error = rep(NA, 2000))
      }
    )
    save(analyse, file = paste0(dir.save, mt, '/', filename))
    
    # load(file = paste0(dir.save, mt, '/', filename))
    
    data.frame(
      method = mt,
      F.set = f,
      boot = b,
      cv = analyse$converge,
      RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
    )
  }
stopCluster(cl)
stopImplicitCluster()
result %>%
  group_by(method, F.set) %>%
  summarise(RMSE = min(RMSE, na.rm = T)) %>%
  mutate(method = factor(method, levels = methods)) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = RMSE, color = method))