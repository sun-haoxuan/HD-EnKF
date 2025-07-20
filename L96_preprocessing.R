rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)

dir.root = '~/L96'
setwd(dir.root)

'./Code/Assimilation' %>% list.files(pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
'./Code/Calculation' %>% list.files(pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)
'./Code/Module' %>% list.files(pattern = "\\.R$", full.names = TRUE) %>% purrr::walk(source)

dir.save.xt = 'Input/TrueState/'
dir.save.yo = 'Input/Observation/'
dir.save.dM = 'Input/distMat/'
dir.save.H = 'Input/ObsMat/'
dir.save.R = 'Input/ObsErr/'

dir.create(dir.save.xt)
dir.create(dir.save.yo)
dir.create(dir.save.dM)
dir.create(dir.save.H)
dir.create(dir.save.R)

### p specific ####
for(p in c(40, 100, 200)){
  option = list(
    boot = NULL,
    ## L96 State
    p = p,
    q = NULL,
    S = 2000,
    S.seq = 4,
    F.true = 8,
    h = 0.05,
    sigma = 1,
    rho = 0.5,
    n = NULL,
    F.set = NULL,
    ## L96 Ensemble
    sigma.in = 0.1,
    method = NULL
  )
  
  x.t = L96_TrueState(option) 
  save(x.t, file = paste0(dir.save.xt, 'L96_p', p, '.Rdata'))
  
  distMat = cal_distMat(option)
  save(distMat, file = paste0(dir.save.dM, 'L96_p', p, '.Rdata'))
}

### p q b specific #####
param1 = param2 = param3 = param4 = param5 = expand.grid(p = c(40, 100, 200), b = 1:2)
param1$q = param1$p
param2$q = param2$p * 0.3
param3$q = param3$p * 0.5
param4$q = param4$p * 0.7
param5$q = param5$p * 0.9
param = param1 %>% rbind(param2) %>% rbind(param3) %>% rbind(param4) %>% rbind(param5)

cores = 10
cl = makeCluster(cores)
registerDoParallel(cl, cores = cores)
foreach(i = 1:nrow(param)) %do%
  {
    p = param$p[i]
    q = param$q[i]
    b = param$b[i]
    
    option = list(
      boot = b,
      ## L96 State
      p = p,
      q = q,
      S.seq = 4,
      S = 2000,
      F.true = 8,
      h = 0.05,
      sigma = 1,
      rho = 0.5,
      n = NULL,
      F.set = NULL
    )
    
    H = cal_H(option)
    save(H, file = paste0(dir.save.H, 'L96_p', p, '_q', q, '_boot', b, '.Rdata'))
    
    Rlist = cal_R(option)
    save(Rlist, file = paste0(dir.save.R, 'L96_p', p, '_q', q, '_boot', b, '.Rdata'))
    
    load(paste0(dir.save.xt, 'L96_p', p, '.Rdata'))
    y.o = L96_Observation(x.t, Rlist$R, option)
    save(y.o, file = paste0(dir.save.yo, 'L96_p', p, '_q', q, '_boot', b, '.Rdata'))
  }
stopCluster(cl)
stopImplicitCluster()
