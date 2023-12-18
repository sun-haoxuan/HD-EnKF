rm(list = ls())
setwd("~/EnKF_GTE")
library(tidyverse)
library(foreach)
library(doParallel)
library(ggsci)

methods = c('standard', 'oracle', 'inflation',  'banding', 'tapering', 'GC')
param = list()
for(mt in methods){
  for(f in 4:12){
    for(b in 1:10){
      param[['f']] = c(param[['f']], f)
      param[['b']] = c(param[['b']], b)
      param[['mt']] = c(param[['mt']], mt)
    }
  }
}
param = as.data.frame(param)

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
  # 
  # rbind(c(100, 100, 20, 4)) %>%
  # rbind(c(100, 100, 30, 4)) %>%
  # rbind(c(100, 100, 40, 4)) %>%
  # 
  # rbind(c(200, 200, 20, 4)) %>%
  # rbind(c(200, 200, 30, 4)) %>%
  # rbind(c(200, 200, 40, 4)) %>%

  slice(-1)

cores = 90
cl = makeCluster(cores)
registerDoParallel(cl, cores = cores)

result.all = foreach(i = 1:nrow(setting), .combine = 'rbind') %do% 
  {
    print(setting[i, ])
    p = setting$p[i]
    q = setting$q[i]
    n = setting$n[i]
    S.seq = setting$S.seq[i]
    dir.save = paste0('Output/L96_param5/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '/')
    
    result = foreach(j = 1:nrow(param), .packages = 'RSpectra', .combine = 'rbind') %do%
      {
        f = param$f[j]
        b = param$b[j]
        mt = param$mt[j]
        filename = paste0(mt, '_F', f, '_B', b, '.Rdata')
        
        load(file = paste0(dir.save, filename))
        
        data.frame(
          method = mt,
          F.set = f,
          boot = b,
          cv = analyse$converge,
          RMSE = sqrt(mean(analyse$error[1001:2000] ^ 2)),
          obj = mean(analyse$likelihood[251:500]),
          iteration = mean(analyse$iteration),
          bandwidth = mean(analyse$bandwidth)
        )
      }
    g = result %>%
      dplyr::select(method, F.set, RMSE) %>%
      mutate(method = factor(method, levels = methods)) %>%
      group_by(method, F.set) %>%
      summarise(m = mean(RMSE, na.rm = T),
                l = min(RMSE, na.rm = T),
                u = max(RMSE, na.rm = T)) %>%
      # arrange(RMSE, by_group = T) %>%
      # slice(1) %>%
      ggplot()+ 
      geom_line(aes(x = F.set, y = m, color = method))+
      geom_ribbon(aes(x = F.set, ymin = l, ymax = u, color = method), alpha = .2)+
      scale_color_d3()+
      labs(x = 'Force term', y = 'RMSE')+
      theme_bw()+
      theme(legend.title = element_blank())
    write.csv(result, paste0('Output/L96_set1_p2/Result/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '.csv'), row.names = F)
    ggsave(paste0('Output/L96_set1_p2/Result/L96_p', p, '_q', q, '_n', n, '_seq', S.seq, '.png'), g, width = 6, height = 4, dpi = 300)
    
    result %>%
      mutate(p = p, q = q, n = n, S.seq = S.seq, .before = 1)
  }
stopCluster(cl)
stopImplicitCluster()
write.csv(result.all, 'Output/L96_set1_p2/Result/Result_all-1120.csv', row.names = F)

result.all %>%
  mutate(method = factor(method, methods)) %>%
  group_by(p, n, F.set, method) %>%
  summarise(m = mean(RMSE, na.rm = T),
            l = min(RMSE, na.rm = T),
            u = max(RMSE, na.rm = T)) %>%
  ggplot()+
  geom_line(aes(x = F.set, y = m, color = method))+
  geom_ribbon(aes(x = F.set, ymin = l, ymax = u, color = method), alpha = .2)+
  facet_wrap(vars(p, n), ncol = 3,
             labeller = label_bquote(p==.(p)~n==.(n)))+
  labs(x = 'Force term', y = 'RMSE')+
  ggsci::scale_color_npg()+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
  