rm(list = ls())
library(tidyverse)
library(foreach)

dir.open = 'Output/Result_SWE_max10/'
result = foreach(f = dir(dir.open), .combine = 'rbind') %do% {
  load(paste0(dir.open, f))
  data.frame(
    method = f,
    RMSE_pre = sqrt(mean(analyse$error[1, 2921:5820] ^ 2)),
    RMSE_post = sqrt(mean(analyse$error[1, 5821:8700] ^ 2))
  )
}

sum_RMSE = function(f){
  load(paste0('Output/Result_SWE_1103/', f))
  tmp = data.frame(
    method = strsplit(f, '_')[[1]][1],
    x = 1:8700,
    # RMSE_all = analyse$error[1, ],
    RMSE_u = analyse$error[2, ],
    RMSE_v = analyse$error[3, ],
    RMSE_h = analyse$error[4, ]
  )
  return(tmp)
}
sum_RMSE('standard_K30000_B1_.Rdata') %>% 
  rbind(sum_RMSE('inflation_K30000_B1_.Rdata')) %>% 
  rbind(sum_RMSE('band-iter_K30000_B1.Rdata')) %>%
  gather(type, value, 3:5) %>%
  ggplot()+
  geom_line(aes(x = x, y = value, color = method))+
  facet_wrap(vars(type), ncol = 1)+
  theme_bw()