rm(list = ls())
setwd("E:/Project/EnKF_GTE")
# setwd("~/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)

source('Code/Functions/cal_HPHR.R')
source('Code/Functions/cal_RMSE.R')
source('Code/Functions/cal_taper_cov_L96.R')
source('Code/Functions/L96_State.R')
source('Code/Functions/L96_Analysis.R')
source('Code/Functions/RK4.R')

methods = c('GC-ev-p')
param = data.frame()
for(mt in methods){
  for(f in 4:12){
    is.iter = T
    for(b in 1:1){
      param = rbind(param, data.frame(
        mt = mt, f = f, b = b, is.iter = is.iter
      ))
    }
  }
}

setting = data.frame(
  p = NA,
  q = NA,
  n = NA,
  S.seq = NA
)
setting = setting %>%
  # rbind(c(40, 40, 20, 4)) %>%
  # rbind(c(40, 40, 30, 4)) %>%
  # rbind(c(40, 40, 40, 4)) %>%
  # rbind(c(100, 100, 20, 4)) %>%
  # rbind(c(100, 100, 30, 4)) %>%
  # rbind(c(100, 100, 40, 4)) %>%
  rbind(c(200, 200, 20, 4)) %>%
  rbind(c(200, 200, 30, 4)) %>%
  rbind(c(200, 200, 40, 4)) %>%
  # rbind(c(500, 500, 20, 4)) %>%
  # rbind(c(500, 500, 30, 4)) %>%
  # rbind(c(500, 500, 40, 4)) %>%
  slice(-1)
task.name = 'L96_set25_GC-ev-p_scaleD'
dir.create(paste0('Output/L96_test/', task.name))

cores = 9
cl = makeCluster(cores)
registerDoParallel(cl, cores = cores)

for(i in 1:nrow(setting)){
  print(setting[i, ])
  p = setting$p[i]
  q = setting$q[i]
  n = setting$n[i]
  S.seq = setting$S.seq[i]
  
  dir.save = paste0('Output/L96_test/', task.name , '/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
  dir.create(dir.save)
  
  result = foreach(j = 1:nrow(param), .packages = 'RSpectra', .combine = 'rbind') %dopar%
    {
      print(param[j, ])
      f = param$f[j]
      b = param$b[j]
      mt = param$mt[j]
      is.iter = param$is.iter[j]
      filename = paste0(mt, '_F', f, '_B', b, '.Rdata')
      
      interval = c((log(p) / n) ^ (-1 / 2) / 5, (log(p) / n) ^ (-1 / 2) * 5)

      set.seed(b)
      option = list(
        p = p,
        S = 2000,
        F.true = 8,
        h = 0.05,
        p.o = sort(sample(p, q)),
        S.o = seq(S.seq, 2000, S.seq),
        sigma = 1,
        rho = 0.5,
        n = n,
        method = mt,
        F.set = f,
        is.iter = is.iter,
        interval = interval,
        sigma.in = 0.1
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
      save(analyse, file = paste0(dir.save, filename))
      
      # load(file = paste0(dir.save, filename))
      
      data.frame(
        method = mt,
        F.set = f,
        boot = b,
        cv = analyse$converge,
        RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2))
      )
    }
  result %>%
    group_by(method, F.set) %>%
    summarise(RMSE = min(RMSE, na.rm = T)) %>%
    mutate(method = factor(method, levels = methods)) %>%
    spread(F.set, RMSE) %>%
    print()
  g = result %>%
    group_by(method, F.set) %>%
    summarise(u = max(RMSE, na.rm = T),
              l = min(RMSE, na.rm = T),
              m = mean(RMSE, na.rm = T)) %>%
    mutate(method = factor(method, levels = methods)) %>%
    ggplot()+
    # geom_ribbon(aes(x = F.set, ymin = l, ymax = u, color = method, fill = method), alpha = .1)+
    geom_line(aes(x = F.set, y = m, color = method))+
    # geom_line(aes(x = F.set, y = l, color = method))+
    labs(main = paste(p, q, n, S.seq))
  print(g)
}
stopCluster(cl)
stopImplicitCluster()
